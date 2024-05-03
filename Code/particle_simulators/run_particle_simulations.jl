#=
    For simulating the particle simulator models. Running the code can take
    quite some time (hence we already have the intermediate results saved :)
    Convenient command-line arguments are:

    `test_all_diffusion` - run all simulations for the particle model with diffusion.
    `test_all_no_diffusion` - run all simulations for the particle model without diffusion.
    `all_illustrations` - create all illustration plots (illustrate simulators for the MS).
=#

using LinearAlgebra
using Plots, StatsPlots
using Random
using StatsBase
using Distributions
using ProgressMeter
using CSV
using DataFrames

mutable struct Particle
    const position::Vector{Float64}
    death_time::Float64
    is_on_membrane::Bool
end

# Convert eucliedan to spherical coordinates
function eu_coordinates_to_spherical(c::Vector{T})::Vector{T} where T<:AbstractFloat
    x, y, z = c
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    ϕ = acos(z / r)
    θ = atan(y, x)
    return [ϕ, θ, r]
end

# From spherical coordinates calculate the corresponding euclidean coordinates.
function spherical_coordinates_to_eu(c::Vector{T})::Vector{T} where T<:AbstractFloat
    ϕ, θ, r = c
    return [r * cos(θ) * sin(ϕ), r * sin(θ) * sin(ϕ), r * cos(ϕ)]
end

# Rescale the vector coord_arr_eu to have length r_new.
function rescale_eu_vector!(c::Vector{T}, r_new::T) where T<:AbstractFloat

    len = norm(c)
    if len == 0
        return
    end
    c .*= r_new / len
end

# On a sphere of radius r compute geodistance between a and b
function compute_geodistance_sphere(a::Vector{T}, b::Vector{T}, R::T)::T where T<:AbstractFloat
    dotprod = dot(a, b) / R^2
    # Numerical noise edge-cases
    if dotprod > 1.0
        dotprod = 1.0
    elseif dotprod < -1.0
        dotprod = -1.0
    end
    return R * acos(dotprod)
end

# dr = distance a point should be pushed
function compute_dr_exo(p::Vector{T}, hit_point::Vector{T},
                        A_exo::T, R::T, alpha::T, γ::T)::T where T<:AbstractFloat

    R_prim = sqrt(R^2 + A_exo/(4*π))
    s = compute_geodistance_sphere(p, hit_point, R)

    arg_arccos = 1.0 - (R/R_prim)^2 * (1.0 - cos(s/R)) - A_exo / (2*π*R_prim^2*alpha)
    # Special case when point is on opposite end of the cell (should then not push)
    if arg_arccos <= -1.0
        return 0.0
    else
        return (1.0 / γ) * (R_prim * acos(arg_arccos) - s)
    end
end

function find_point_exocytosis(hit_point::Vector{T},
                               b_point::Vector{T},
                               A_exo::T,
                               R::T; α=0.5, γ=1.0)::Vector{T} where T<:AbstractFloat



    # Check special cases with the opposite end of the cell
    d_bc = compute_dr_exo(hit_point, b_point, A_exo, R, α, γ)
    if d_bc == 0.0
        return deepcopy(b_point)
    end

    a1, a2, a3 = hit_point
    b1, b2, b3 = b_point

    # Compute distances between points
    d_ab = compute_geodistance_sphere(hit_point, b_point, R)
    d_ac = d_bc + d_ab
    ϕ_ac = d_ac / R
    ϕ_bc = d_bc / R

    # All following computations are derived from symbolic solver
    A1 = R^2 * cos(ϕ_ac)
    A2 = R^2 * cos(ϕ_bc)

    nominator = (a1*b2 - a2*b1)
    c1_tmp1 = (A1*b2 - A2*a2) / nominator
    c1_tmp2 = (a2*b3 - a3*b2)
    c1_tmp3 = nominator
    c2_tmp1 = (-A1*b1 + A2*a1) / nominator
    c2_tmp2 = (-a1*b3 + a3*b1)
    c2_tmp3 = nominator

    c3 = c1_tmp3*c2_tmp3*(-c1_tmp1*c1_tmp2*c2_tmp3 - c1_tmp3*c2_tmp1*c2_tmp2) / (c1_tmp2^2*c2_tmp3^2 + c1_tmp3^2*c2_tmp2^2 + c1_tmp3^2*c2_tmp3^2)
    if nominator == 0.0
        c1, c2 = 0.0, 0.0
    else
        c1 = (A1*b2 - A2*a2 + a2*b3*c3 - a3*b2*c3)/(a1*b2 - a2*b1)
        c2 = (-A1*b1 + A2*a1 - a1*b3*c3 + a3*b1*c3)/(a1*b2 - a2*b1)
    end

    return [c1, c2, c3]
end

function get_mean_distance(points, R)
    all_dist = Vector{Float64}(undef, size(points)[2])
    for i in 1:size(points)[2]-1
        p_ref = points[:, i]
        compare_to = points[:, i+1:end]
        try
            all_dist[i] = minimum(R .* acos.(compare_to' * p_ref ./ R^2))
        catch
            all_dist[i] = minimum([compute_geodistance_sphere(p_ref, compare_to[:, k], R) for k in size(compare_to)[2]])
        end
    end
    return mean(all_dist)
end

function find_point_further_away(p_centre, point_furthest_away, R)

    _p_centre = normalize(p_centre)
    _point_furthest_away = normalize(point_furthest_away)

    θ = dot(_p_centre, _point_furthest_away)
    ϕ = ifelse(θ > 0.0, 0.02, -0.02)
    n = cross(_p_centre, _point_furthest_away)

    # Following equations are dervied from a symbolic solver
    A1 = cos(θ + ϕ)
    A2 = 0
    b1, b2, b3 = _point_furthest_away
    a1, a2, a3 = _p_centre
    n1, n2, n3 = n

    c3_1 = (-A1*(a1*n1*n3 + a2*n2*n3 - a3*n1^2 - a3*n2^2) + (-a1*n2 + a2*n1)*sqrt(-A1^2*n1^2 - A1^2*n2^2 - A1^2*n3^2 + a1^2*n2^2 + a1^2*n3^2 - 2*a1*a2*n1*n2 - 2*a1*a3*n1*n3 + a2^2*n1^2 + a2^2*n3^2 - 2*a2*a3*n2*n3 + a3^2*n1^2 + a3^2*n2^2))/(a1^2*n2^2 + a1^2*n3^2 - 2*a1*a2*n1*n2 - 2*a1*a3*n1*n3 + a2^2*n1^2 + a2^2*n3^2 - 2*a2*a3*n2*n3 + a3^2*n1^2 + a3^2*n2^2)
    c3_2 = (-A1*(a1*n1*n3 + a2*n2*n3 - a3*n1^2 - a3*n2^2) + (a1*n2 - a2*n1)*sqrt(-A1^2*n1^2 - A1^2*n2^2 - A1^2*n3^2 + a1^2*n2^2 + a1^2*n3^2 - 2*a1*a2*n1*n2 - 2*a1*a3*n1*n3 + a2^2*n1^2 + a2^2*n3^2 - 2*a2*a3*n2*n3 + a3^2*n1^2 + a3^2*n2^2))/(a1^2*n2^2 + a1^2*n3^2 - 2*a1*a2*n1*n2 - 2*a1*a3*n1*n3 + a2^2*n1^2 + a2^2*n3^2 - 2*a2*a3*n2*n3 + a3^2*n1^2 + a3^2*n2^2)

    c1_1 = (A1*n2 + a2*c3_1*n3 - a3*c3_1*n2)/(a1*n2 - a2*n1)
    c1_2 = (A1*n2 + a2*c3_2*n3 - a3*c3_2*n2)/(a1*n2 - a2*n1)

    c2_1 = (-A1*n1 - a1*c3_1*n3 + a3*c3_1*n1)/(a1*n2 - a2*n1)
    c2_2 =(-A1*n1 - a1*c3_2*n3 + a3*c3_2*n1)/(a1*n2 - a2*n1)

    c_1 = [c1_1, c2_1, c3_1] .* R
    c_2 = [c1_2, c2_2, c3_2] .* R

    dist1 = compute_geodistance_sphere(p_centre, c_1, R)
    dist2 = compute_geodistance_sphere(p_centre, c_2, R)

    if dist1 < dist2
        return c_1
    else
        return c_2
    end
end

function compute_ring_size(points, R; print::Bool=false)

    # I need the ring counting algorithm now
    i_use = findall(x -> x == true, [sum(isnan.(points[:, i])) == 0 for i in 1:size(points)[2]])
    _points = points[:, i_use]
    p_centre = normalize(mean(_points, dims=2)[:]) * R

    distance_to_centre = [compute_geodistance_sphere(p_centre, _points[:, i], R) for i in 1:length(i_use)]
    point_furthest_away = points[:, argmax(distance_to_centre)]
    point_further_away = find_point_further_away(p_centre, point_furthest_away, R)

    θ_list = collect(0:0.02:2*π-0.01)
    v_rot = Vector{Vector{Float64}}(undef, length(θ_list))
    for (i, θ) in pairs(θ_list)
        __point_further_away = normalize(point_further_away)
        __point_centre =  normalize(p_centre)
        v_rot[i] = __point_further_away .* cos(θ) + cross(__point_centre, __point_further_away).*sin(θ) + __point_centre .* dot(__point_centre, __point_further_away) .* (1 - cos(θ))
    end
    for i in eachindex(v_rot)
        v_rot[i] .*= R
    end

    mean_dist = get_mean_distance(_points, R)
    is_inner_point = Vector{Bool}(undef, size(_points)[2]) .* false
    is_outer_point = Vector{Bool}(undef, size(_points)[2]) .* false
    t_values = collect(1e-4:1e-3:1-1e-4)

    local dist_to_all_points
    for i in eachindex(v_rot)
        p_end = v_rot[i]
        for (j, t) in pairs(t_values)
            R_point = normalize(p_centre + t.*(p_end - p_centre)) .* R

            try
                dist_to_all_points = R .* acos.(_points' * R_point ./ R^2)
            catch
                dist_to_all_points = [compute_geodistance_sphere(R_point, _points[:, k], R) for k in size(_points)[2]]
            end

            if minimum(dist_to_all_points) < mean_dist*2.0
                is_inner_point[argmin(dist_to_all_points)] = true
                break
            end
        end
        for (j, t) in pairs(reverse(t_values))
            R_point = normalize(p_centre + t.*(p_end - p_centre)) .* R

            try
                dist_to_all_points = R .* acos.(_points' * R_point ./ R^2)
            catch
                dist_to_all_points = [compute_geodistance_sphere(R_point, _points[:, k], R) for k in size(_points)[2]]
            end

            if minimum(dist_to_all_points) < mean_dist*2.0
                is_outer_point[argmin(dist_to_all_points)] = true
                break
            end
        end
    end

    inner_points = findall(x -> x == true, is_inner_point)
    outer_points = findall(x -> x == true, is_outer_point)
    dist_inner_points = [compute_geodistance_sphere(p_centre, _points[:, k], R) for k in inner_points]
    dist_outer_points = [compute_geodistance_sphere(p_centre, _points[:, k], R) for k in outer_points]

    colors = Vector{String}(undef, length(inner_points))
    for i in eachindex(colors)
        if is_inner_point[i] == true
            colors[i] = "blue"
        elseif is_outer_point[i] == true
            colors[i] = "orange"
        else
            colors[i] = "white"
        end
    end
    p = scatter3d(_points[1, inner_points], _points[2, inner_points], _points[3, inner_points], color="blue",
              xlim=[-R, R].*1.2,
              ylim=[-R, R].*1.2,
              zlim=[-R, R].*1.2,
              xlab="x", ylab="y", zlab="z",
              camera=(45, 10))
    p = scatter3d!(_points[1, outer_points], _points[2, outer_points], _points[3, outer_points], color="orange")
    if print == true
        display(p)
    end

    return dist_inner_points, dist_outer_points, p
end

# I want to distribute points around another point
# frac_cover - fraction of membrane around the point allowed to cover
function generate_points_around_another!(points::Matrix{T},
                                         point_centre::Vector{T},
                                         frac_cover::T,
                                         R::T) where T<:AbstractFloat

    # Generating point
    d_cutoff = R * acos(1 - 2 * frac_cover)
    n_points = size(points)[2]
    k = 1
    while k ≤ n_points
        for i in frac_cover * 100 * 1000
            point = normalize(rand(Uniform(-1, 1), 3)) * R
            dist = compute_geodistance_sphere(point, point_centre, R)
            if dist < d_cutoff
                points[:, k] .= point
                k += 1
                break
            end
        end
    end
end

function simulate_model(R, point_centre, A_exo, n_steps,
                        power_weight, frac_cover; print::Bool=false,
                        cables::Bool=false,
                        deliver_particles::Bool=false,
                        n_points::Int64=1000,
                        α=0.5, γ=1.0)

    points = Matrix{Float64}(undef, (3, 1000))

    # Generating point
    generate_points_around_another!(points, point_centre, frac_cover, R)

    p = plot = scatter3d(points[1, :], points[2, :], points[3, :],
            xlim=[-R, R].*1.2,
            ylim=[-R, R].*1.2,
            zlim=[-R, R].*1.2,
            xlab="x", ylab="y", zlab="z",
            camera=(45, 10))
    p = plot = scatter3d!([point_centre[1]], [point_centre[2]], [point_centre[3]],
            color="red")
    if print == true
        display(p)
    end

    points_start = deepcopy(points)
    # Get distance to center
    dist_to_centre = [compute_geodistance_sphere(point_centre, points_start[:, i], R) for i in 1:n_points]
    weights = exp10.(log10(1) .- power_weight .* log10.(dist_to_centre .+ 1e-8))
    weights = ProbabilityWeights( weights ./ sum(weights))

    local p
    local cables_i
    for step in 1:n_steps

        if (step == 1 || step % 10 == 0) && cables == true
            cables_i = sample(1:n_points, weights, 10)
        end

        if cables == false
            new_point_centre = points_start[:, sample(1:n_points, weights)]
        else
            i_centre = sample(cables_i)
            new_point_centre = points_start[:, i_centre]
        end

        # If new particles should be inserted
        if deliver_particles == true
            delivered_particles = Matrix{Float64}(undef, 3, 10)
            generate_points_around_another!(delivered_particles, new_point_centre, frac_cover / 10, R)
            points = hcat(points, delivered_particles)
        end

        new_points = similar(points)
        for i in 1:size(points)[2]
            new_points[:, i] = find_point_exocytosis(new_point_centre, points[:, i], A_exo, R,
                                                     α=α, γ=γ)
        end
        points = new_points

        p = scatter3d(points[1, :], points[2, :], points[3, :],
                xlim=[-R, R].*1.2,
                ylim=[-R, R].*1.2,
                zlim=[-R, R].*1.2,
                xlab="x", ylab="y", zlab="z",
                camera=(45, 10))

        p = scatter3d!([new_point_centre[1]], [new_point_centre[2]], [new_point_centre[3]],
            color="red")

        if print == true
            display(p)
            sleep(0.2)
        end
    end

    return p, points
end

function run_simulation(focus_center_list, n_repeat::Int64, fraction_start,
                        iterations_simulate, name_save::String;
                        cables::Bool=false,
                        deliver_particles::Bool=false,
                        α::Float64=0.5,
                        γ::Float64=1.0)

    R = 2.5
    point_centre = normalize([1.0, -1.0, 0.0]) .* R
    r_exo = 0.05
    A_exo = 4π*r_exo^2

    # Run point simulations
    _focus_center_list = repeat(focus_center_list, inner=n_repeat)
    points_list = Vector{Matrix{Float64}}(undef, length(_focus_center_list))
    points_figures = Vector{Any}(undef, length(_focus_center_list))
    @showprogress 1 "Simulating points..." for i in eachindex(_focus_center_list)
        points_figures[i], points_list[i] = simulate_model(R, point_centre, A_exo, iterations_simulate, _focus_center_list[i],
                                                           fraction_start, cables=cables, deliver_particles=deliver_particles,
                                                           α=α, γ=γ)
    end

    # Process densities
    inner_points = Vector{Vector{Float64}}(undef, length(_focus_center_list))
    outer_points = Vector{Vector{Float64}}(undef, length(_focus_center_list))
    inner_outer_figures = Vector{Any}(undef, length(_focus_center_list))
    @showprogress 1 "Computing ring sizes..." for i in eachindex(points_list)
        inner_points[i], outer_points[i], inner_outer_figures[i] = compute_ring_size(points_list[i], R)
    end

    # Export Data to CSV file (and do the plotting in R)
    df = DataFrame()
    for i in eachindex(inner_points)

        inner, outer = inner_points[i], outer_points[i]
        # A ring failed to form. If particles are delivered inner diameter is
        # not checked
        if (length(inner) < 20 && deliver_particles == false) || length(outer) < 20
            _df1 = DataFrame(id = i,
                            strength = _focus_center_list[i],
                            distance = NaN,
                            point_type = "inner")
            _df2 = DataFrame(id = i,
                            strength = _focus_center_list[i],
                            distance = NaN,
                            point_type = "outer")
            append!(df, _df1)
            append!(df, _df2)
            continue
        end
        if deliver_particles == false
            _df1 = DataFrame(id = i,
                            strength = _focus_center_list[i],
                            distance = inner,
                            point_type = "inner")
        else
            _df1 = DataFrame()
        end
        _df2 = DataFrame(id = i,
                        strength = _focus_center_list[i],
                        distance = outer,
                        point_type = "outer")
        append!(df, _df1)
        append!(df, _df2)
    end

    dir_save = joinpath(@__DIR__, "..", "..", "Intermediate", "Simulations", "Particle_no_diff")
    if !isdir(dir_save)
        mkpath(dir_save)
    end
    file_path = joinpath(dir_save, name_save * ".csv")
    @info "Saving results at $file_path"
    CSV.write(file_path, df)
end

function update_pos!(particle::AbstractVector{Float64}, diffusion_strength::Float64, R::Float64)
    # diffusion_strength sqrt(2*D*Δt)
    particle .+= diffusion_strength .* randn(3)
    normalize!(particle)
    particle .*= R
end

# Can be modelled with a classical Gillespie algorithm
function compute_membrane_time(k_off::Float64, t::Float64)::Float64
    return -log(rand()) / k_off + t
end

function deliver_new_particles(particles::Vector{Particle},
                               n_particles_deliver::Int64,
                               point_centre::Vector{Float64},
                               t::Float64,
                               frac_cover::Float64,
                               k_off::Float64,
                               R::Float64)::Vector{Particle}

    delivered_particles = Matrix{Float64}(undef, 3, n_particles_deliver)
    generate_points_around_another!(delivered_particles, point_centre, frac_cover / 10, R)
    n_points_assigned = 0
    for ix in eachindex(particles)
        if particles[ix].is_on_membrane == false
            particles[ix].is_on_membrane = true
            particles[ix].position .= delivered_particles[:, n_points_assigned+1]
            particles[ix].death_time = compute_membrane_time(k_off, t)
            n_points_assigned += 1
        end
        if n_points_assigned == n_particles_deliver
            break
        end
    end

    new_list = [Particle(Vector{Float64}(undef, 3), Inf, false) for i in 1:(n_particles_deliver-n_points_assigned)]
    for ix in eachindex(new_list)
        new_list[ix].is_on_membrane = true
        new_list[ix].position .= delivered_particles[:, n_points_assigned+1]
        new_list[ix].death_time = compute_membrane_time(k_off, t)
        n_points_assigned += 1
        if n_points_assigned == n_particles_deliver
            break
        end
    end

    return vcat(particles, new_list)
end

function recruit_particles(Δt::T, k_recruit::T, R::T, frac_cover::T, k_off::T, t::T, point_centre::Vector{T})::Vector{Particle} where T<:AbstractFloat

    n_recruited = rand(Poisson(k_recruit*Δt))

    # Randomly select location within what is covered
    new_coord = Matrix{Float64}(undef, 3, n_recruited)
    generate_points_around_another!(new_coord, point_centre, frac_cover, R)

    new_particles = [Particle(Vector{Float64}(undef, 3), Inf, false) for i in 1:n_recruited]
    for ix in eachindex(new_particles)
        new_particles[ix].is_on_membrane = true
        new_particles[ix].position .= new_coord[:, ix]
        new_particles[ix].death_time = compute_membrane_time(k_off, t)
    end

    return new_particles
end

function simulate_particles(Δt::Float64,
                            frac_cover::Float64,
                            n_particles_deliver::Int64,
                            k_off::Float64,
                            power_weight::Float64,
                            t_end::Float64,
                            λhit::Float64,
                            n_save::Int64,
                            Dm::Float64,
                            R::Float64,
                            α::Float64,
                            γ::Float64,
                            k_recruit::Float64,
                            have_exo::Bool,
                            recruit_with_exo::Bool,
                            recruit_particles_reaction::Bool,
                            use_cables::Bool)

    time_hit::Float64 = -log(rand()) * λhit

    data_save = Vector{Matrix{Float64}}(undef, n_save)
    n_steps = length(0.0:Δt:t_end)
    i_save = ceil.(collect(1:n_steps ./ n_save:n_steps)) |> Vector{Int64}

    point_centre = normalize([1.0, -1.0, 0.0]) .* R
    n_points = 10
    points = Matrix{Float64}(undef, 3, n_points)
    generate_points_around_another!(points, point_centre, frac_cover, R)
    particles = Vector{Particle}(undef, n_points)
    for i in eachindex(particles)
        particles[i] = Particle(points[:, i], compute_membrane_time(k_off, 0.0), true)
    end

    # For generating delivery points
    points_start = deepcopy(points)
    dist_to_centre = [compute_geodistance_sphere(point_centre, points_start[:, i], R) for i in 1:n_points]
    weights = exp10.(log10(1) .- power_weight .* log10.(dist_to_centre .+ 1e-8))
    weights = ProbabilityWeights( weights ./ sum(weights))

    λ_cable = 1 / 60.0
    if use_cables == true
        _cables_i = sample(1:n_points, weights, 10)
        cables = Vector{Vector{Float64}}(undef, 10)
        for i in eachindex(cables)
            cables[i] = [_cables_i[i], -log(rand()) / λ_cable]
        end
    end

    # Pre-compute what can be pre-computed
    diffusion_strength = sqrt(2*Dm*Δt)
    t::Float64 = 0.0
    j = 1
    for it in 1:n_steps

        # Check if should update cables
        if use_cables == true
            for i in eachindex(cables)
                if cables[i][2] > t
                    cables[i][2] = -log(rand()) / λ_cable + t
                end
            end
        end

        # Exocytosis occured!
        if time_hit < t

            if use_cables == false
                new_point_centre = points_start[:, sample(1:n_points, weights)]
            else
                i_centre = Int64(sample([cables[k][1] for k in eachindex(cables)]))
                new_point_centre = points_start[:, i_centre]
            end

            if have_exo == true
                for ix in eachindex(particles)
                    particles[ix].position .= find_point_exocytosis(new_point_centre, particles[ix].position, A_exo, R,
                                                                    α=α, γ=γ)
                end
            end
            time_hit = t - log(rand())*λhit

            if recruit_with_exo == true
                particles = deliver_new_particles(particles, n_particles_deliver, new_point_centre, t, frac_cover, k_off, R)
            end
        end

        if recruit_particles_reaction == true
            __point_centre = normalize([1.0, -1.0, 0.0]) .* R
            new_particles = recruit_particles(Δt, k_recruit, R, frac_cover, k_off, t, __point_centre)
            particles = vcat(particles, new_particles)
        end

        # Check if any membrane bound molecules are out of comission
        for ix in eachindex(particles)
            if particles[ix].is_on_membrane == false
                continue
            end
            if particles[ix].death_time < t
                particles[ix].is_on_membrane = false
            end
        end

        # Expansive allocation, should not happen to often
        if it ∈ i_save
            n_particles_membrane::Int64 = 0
            for ix in eachindex(particles)
                n_particles_membrane += particles[ix].is_on_membrane
            end
            _matrix_save = Matrix{Float64}(undef, 3, n_particles_membrane)
            k = 1
            for ix in eachindex(particles)
                if particles[ix].is_on_membrane == true
                    _matrix_save[:, k] .= particles[ix].position
                    k += 1
                end
            end
            data_save[j] = _matrix_save
            j += 1
        end

        for ix in eachindex(particles)
            if particles[ix].is_on_membrane == true
                update_pos!(particles[ix].position, diffusion_strength, R)
            end
        end

        t += Δt
    end
    return data_save
end

function compute_dist_to_centre(point_centre::Vector{T}, points::Matrix{T}, R::T)::Vector{T} where T<:AbstractFloat
    dist_to_centre = Vector{Float64}(undef, size(points)[2])
    for i in eachindex(dist_to_centre)
        dist_to_centre[i] = compute_geodistance_sphere(points[:, i], point_centre, R)
    end
    return dist_to_centre
end

function run_simulate_particles(Δt::Float64,
                                tag_save::String,
                                n_repeat::Int64,
                                frac_cover::Float64,
                                n_particles_deliver::Int64,
                                k_off::Float64,
                                power_weight_list::Vector{Float64};
                                t_end::Float64=400.0,
                                λhit::Float64=0.4,
                                n_save::Int64=50,
                                Dm::Float64=0.01,
                                R::Float64=2.5,
                                α::Float64=0.5,
                                γ::Float64=1.0,
                                k_recruit::Float64=40.0,
                                have_exo::Bool=true,
                                recruit_with_exo::Bool=true,
                                recruit_particles_reaction::Bool=false,
                                use_cables::Bool=false)


    power_weight_list = repeat(power_weight_list, n_repeat)
    data_save = Vector{Vector{Matrix{Float64}}}(undef, length(power_weight_list))
    @showprogress 1 "Simulating system..." for i in eachindex(data_save)
        data_save[i] = simulate_particles(Δt, frac_cover, n_particles_deliver, k_off, power_weight_list[i], t_end, λhit,
                                          n_save, Dm, R, α, γ, k_recruit, have_exo, recruit_with_exo, recruit_particles_reaction, use_cables)
    end

    point_centre = normalize([1.0, -1.0, 0.0]) .* R
    dist_to_centre = Vector{Vector{Float64}}(undef, length(data_save))
    for i in eachindex(data_save)
        dist_to_centre[i] = compute_dist_to_centre(point_centre, data_save[i][end], R)
    end

    df = DataFrame()
    for i in eachindex(dist_to_centre)
        _df = DataFrame(runindex = i,
                        dist = dist_to_centre[i],
                        weight_centre = power_weight_list[i])
        append!(df, _df)
    end

    dir_save = joinpath(@__DIR__, "..", "..", "Intermediate", "Simulations", "Particle_with_diff")
    if recruit_with_exo == true && use_cables == false
        dir_save = joinpath(dir_save, "With_diffusion")
    elseif recruit_with_exo == true && use_cables == true
        dir_save = joinpath(dir_save, "With_diffusion_cables")
    else
        dir_save = joinpath(dir_save, "With_diffusion_reac_rec")
    end
    if !isdir(dir_save)
        mkpath(dir_save)
    end
    path_save = joinpath(dir_save, tag_save * ".csv")
    CSV.write(path_save, df, append=isfile(path_save))

    return dist_to_centre
end

# Make images to plot the results
function sphere(r, C)   # r: radius; C: center [cx,cy,cz]
    n = 100
    u = range(-π, π; length = n)
    v = range(0, π; length = n)
    x = C[1] .+ r*cos.(u) * sin.(v)'
    y = C[2] .+ r*sin.(u) * sin.(v)'
    z = C[3] .+ r*ones(n) * cos.(v)'
    return x, y, z
end

# General parameters
R = 2.5
point_centre = normalize([1.0, -1.0, 0.0]) .* R
r_exo = 0.05
A_exo = 4π*r_exo^2

if ARGS[1] == "Dm_particles" || ARGS[1] == "test_all_diffusion"
    # Dm experiment
    k_off = 0.1
    Δt = 0.01

    Dm = parse(Float64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])
    use_cables = parse(Bool, ARGS[4])

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "Dm" * string(Dm) * "frac_cover" * string(frac_cover)
    if use_cables == true
        _tag_save *= "cables"
    end
    run_simulate_particles(Δt, _tag_save, 10, frac_cover, 100, k_off, centred_list; Dm=Dm, use_cables=use_cables)
end

if ARGS[1] == "Dm_particles_no_exo" || ARGS[1] == "test_all_diffusion"
    # Dm experiment
    k_off = 0.1
    Δt = 0.01

    Dm = parse(Float64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])
    use_cables = parse(Bool, ARGS[4])

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "Dm" * string(Dm) * "noexo_frac_cover" * string(frac_cover)
    if use_cables == true
        _tag_save *= "cables"
    end
    run_simulate_particles(Δt, _tag_save, 10, frac_cover, 100, k_off, centred_list; Dm=Dm, have_exo=false, use_cables=use_cables)
end

if ARGS[1] == "Amount_delivered" || ARGS[1] == "test_all_diffusion"
    # Dm experiment
    k_off = 0.1
    Δt = 0.01
    Dm = 0.0045
    amount = parse(Int64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])
    use_cables = parse(Bool, ARGS[4])

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "ammount" * string(amount) * "frac_cover" * string(frac_cover)
    if use_cables == true
        _tag_save *= "cables"
    end
    run_simulate_particles(Δt, _tag_save, 10, frac_cover, amount, k_off, centred_list; Dm=Dm, use_cables=use_cables)
end

if ARGS[1] == "Koff" || ARGS[1] == "test_all_diffusion"
    # Dm experiment
    Δt = 0.01
    Dm = 0.0045
    k_off = parse(Float64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])
    use_cables = parse(Bool, ARGS[4])

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "koff" * string(k_off) * "frac_cover" * string(frac_cover)
    if use_cables == true
        _tag_save *= "cables"
    end
    run_simulate_particles(Δt, _tag_save, 10, frac_cover, 100, k_off, centred_list; Dm=Dm, use_cables=use_cables)
end

if ARGS[1] == "Koff_low_diffusion" || ARGS[1] == "test_all_diffusion"
    # Dm experiment
    Δt = 0.01
    Dm = 0.00045
    k_off = parse(Float64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])
    use_cables = parse(Bool, ARGS[4])

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "koff" * string(k_off) * "low_diffusionDm" * string(Dm) * "frac_cover" * string(frac_cover)
    if use_cables == true
        _tag_save *= "cables"
    end
    run_simulate_particles(Δt, _tag_save, 10, frac_cover, 100, k_off, centred_list; Dm=Dm, use_cables=use_cables)
end

if ARGS[1] == "Koff_low_diffusion_no_exo" || ARGS[1] == "test_all_diffusion"
    # Dm experiment
    Δt = 0.01
    Dm = 0.00045
    k_off = parse(Float64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])
    use_cables = parse(Bool, ARGS[4])

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "koff" * string(k_off) * "low_diffusionDm" * string(Dm) * "frac_cover" * string(frac_cover)
    if use_cables == true
        _tag_save *= "cables"
    end
    run_simulate_particles(Δt, _tag_save, 10, frac_cover, 100, k_off, centred_list; Dm=Dm, have_exo=false, use_cables=use_cables)
end

if ARGS[1] == "Amount_delivered_low_diffusion" || ARGS[1] == "test_all_diffusion"
    # Dm experiment
    k_off = 0.1
    Δt = 0.01
    Dm = 0.00045
    amount = parse(Int64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])
    use_cables = parse(Bool, ARGS[4])

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "ammount" * string(amount) * "low_diffusionDm" * string(Dm) * "frac_cover" * string(frac_cover)
    if use_cables == true
        _tag_save *= "cables"
    end
    run_simulate_particles(Δt, _tag_save, 10, frac_cover, amount, k_off, centred_list; Dm=Dm, use_cables=use_cables)
end

if ARGS[1] == "Test_size_pole" || ARGS[1] == "test_all_no_diffusion"

    frac_cover_list = [0.01, 0.02, 0.03, 0.04, 0.05]
    α_list = [0.4, 0.5, 0.8]

    for frac_cover in frac_cover_list
        for α in α_list
            @info "Testing with frac_cover = $frac_cover and α = $α"
            name_save = "Frac_cover" * string(frac_cover) * "_alpha" * string(α)
            run_simulation(collect(0.3:0.1:5.0), 5, frac_cover, 40, name_save, α=α)
        end
    end
end

if ARGS[1] == "Diffusion_recruitment_Dm" || ARGS[1] == "test_all_diffusion"

    # Dm experiment
    k_off = 0.1
    Δt = 0.01

    Dm = parse(Float64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "Dm" * string(Dm) * "frac_cover" * string(frac_cover)
    run_simulate_particles(Δt, _tag_save, 10, frac_cover, 100, k_off, centred_list;
                           Dm=Dm, recruit_with_exo=false, recruit_particles_reaction=true, k_recruit=50.0, t_end=800.0)
end

if ARGS[1] == "Diffusion_recruitment_alpha" || ARGS[1] == "test_all_diffusion"

    @info "Running Diffusion_recruitment_alpha"

    # Dm experiment
    k_off = 0.1
    Δt = 0.01
    Dm = 0.00045

    α = parse(Float64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "alpha" * string(α) * "frac_cover" * string(frac_cover)
    run_simulate_particles(Δt, _tag_save, 1, frac_cover, 100, k_off, centred_list;
                           Dm=Dm, recruit_with_exo=false, recruit_particles_reaction=true, k_recruit=50.0, α=α, t_end=800.0)
end

if ARGS[1] == "Diffusion_recruitment_koff" || ARGS[1] == "test_all_diffusion"

    @info "Running Diffusion_recruitment_koff"

    # Dm experiment
    Δt = 0.01
    Dm = 0.00045

    k_off= parse(Float64, ARGS[2])
    frac_cover = parse(Float64, ARGS[3])
    α = 0.5

    centred_list = collect(0.3:0.2:3.0)
    _tag_save = "koff" * string(k_off) * "frac_cover" * string(frac_cover)
    run_simulate_particles(Δt, _tag_save, 1, frac_cover, 100, k_off, centred_list;
                           Dm=Dm, recruit_with_exo=false, recruit_particles_reaction=true, k_recruit=50.0, α=α, t_end=800.0)
end

if ARGS[1] == "Test_cables" || ARGS[1] == "test_all_no_diffusion"

    frac_cover_list = [0.01, 0.02, 0.03, 0.04, 0.05]
    α_list = [0.4, 0.5, 0.8]

    for frac_cover in frac_cover_list
        for α in α_list
            @info "Testing cables with frac_cover = $frac_cover and α = $α"
            name_save = "Cables_frac_cover" * string(frac_cover) * "_alpha" * string(α)
            run_simulation(collect(0.3:0.1:5.0), 5, frac_cover, 40, name_save, α=α, cables=true)
        end
    end
end

if ARGS[1] == "Test_delivery" || ARGS[1] == "test_all_no_diffusion"

    frac_cover_list = [0.01, 0.02, 0.03, 0.04, 0.05]
    α_list = [0.4, 0.5, 0.8]

    for frac_cover in frac_cover_list
        for α in α_list
            @info "Testing delivery with frac_cover = $frac_cover and α = $α"
            name_save = "Deliver_frac_cover" * string(frac_cover) * "_alpha" * string(α)
            run_simulation(collect(0.3:0.1:5.0), 8, frac_cover, 40, name_save, α=α, deliver_particles=true)
        end
    end
end

if ARGS[1] == "Illustration_particle_diffusion" || ARGS[1] == "all_illustrations"
    #=
        Simulate particle model with diffusion
    =#
    # Parameter needed by simulator
    dir_save = joinpath(@__DIR__, "..", "..", "Results", "Particle_simulator_diff", "Illustration")
    if !isdir(dir_save)
        mkpath(dir_save)
    end
    Δt = 0.01
    frac_cover=0.01
    n_particles_deliver=100
    k_off = 1.0
    power_weight = 3.0
    t_end = 400.0
    λhit=0.4
    n_save=25
    R=2.5
    α=0.2
    γ=1.0
    have_exo=true
    recruit_with_exo = false
    recruit_particles_reaction = true

    k_off = 0.1
    k_recruit = 60.0
    data_save = Vector{Vector{Matrix{Float64}}}(undef, 2)

    # Fast Dm
    Dm = 0.00005
    Random.seed!(123)
    res = simulate_particles(Δt, frac_cover, n_particles_deliver, k_off, power_weight,
                            t_end, λhit, n_save, Dm, R, α, γ, k_recruit, have_exo, recruit_with_exo, recruit_particles_reaction, false)

    i_list = [1, 2, 4, 5, 10]
    for i in i_list
        _t = (t_end - 0.0) / n_save * (i - 1)
        dist_to_centre = compute_dist_to_centre(point_centre, res[i], R)
        p1 = surface(sphere(2.5, [0.0, 0.0, 0.0]),
                    camera=(45, 10),
                    xlim=[-R, R].*1.2,
                    ylim=[-R, R].*1.2,
                    zlim=[-R, R].*1.2,
                    fill="grey",
                    xlab="x", ylab="y", zlab="z",
                    title = "t = $_t [s]",
                    legend=false)
        p1 = scatter3d!(res[i][1, :], res[i][2, :], res[i][3, :],
                        alpha=0.7,
                        camera=(45, 10))
        path_save = joinpath(dir_save, "t$_t.svg")
        savefig(p1, path_save)

        p2 = density(dist_to_centre, title = "t = $_t", label=false)
        path_save = joinpath(dir_save, "Dens_t$_t.svg")
        savefig(p2, path_save)
    end
end

if ARGS[1] == "Illustration_particle_no_diffusion" || ARGS[1] == "all_illustrations"
    # Parameter needed by simulator
    dir_save = joinpath(@__DIR__, "..", "..", "Results", "Particle_simulator_no_diff", "Illustration")
    if !isdir(dir_save)
        mkpath(dir_save)
    end

    power_weight = 3.0
    frac_cover = 0.01
    α=0.5
    γ=1.0
    _ilist = [1, 5, 10, 15, 20]
    for i in _ilist
        p, points = simulate_model(R, point_centre, A_exo, i,
                                  power_weight, frac_cover; print=false,
                                  cables=false,
                                  deliver_particles=false,
                                  n_points=100,
                                  α=α, γ=γ)

        p1 = surface(sphere(2.5, [0.0, 0.0, 0.0]),
                    camera=(45, 10),
                    xlim=[-R, R].*1.2,
                    ylim=[-R, R].*1.2,
                    zlim=[-R, R].*1.2,
                    fill="grey",
                    xlab="x", ylab="y", zlab="z",
                    title = "Iteration = $i",
                    legend=false)
        p1 = scatter3d!(points[1, :], points[2, :], points[3, :],
                        alpha=0.7,
                        camera=(45, 10))
        path_save = joinpath(dir_save, "it$i.svg")
        savefig(p1, path_save)
    end
end

#=
    Plotting the sampling if the initial values is likelly best done via a density plot I think,
    then we can work from that point
=#
if ARGS[1] == "Data_illustration_hit" || ARGS[1] == "all_illustrations"

    dir_save = joinpath(@__DIR__, "..", "..", "Results", "Particle_simulator_no_diff", "Hit_dist")
    if !isdir(dir_save)
        mkpath(dir_save)
    end

    Random.seed!(123)
    n_points = 100000
    weight_list = [0.1, 0.5, 1.0, 1.5, 2.0]
    for power_weight in weight_list
        n_points = 1000
        frac_cover = 0.01
        points = Matrix{Float64}(undef, (3, n_points))
        # Generating point
        generate_points_around_another!(points, point_centre, frac_cover, R)

        points_start = deepcopy(points)
        # Get distance to center
        dist_to_centre = [compute_geodistance_sphere(point_centre, points_start[:, i], R) for i in 1:n_points]
        weights = exp10.(log10(1) .- power_weight .* log10.(dist_to_centre .+ 1e-8))
        weights = ProbabilityWeights( weights ./ sum(weights))

        i_sample = sample(1:n_points, weights, n_points; replace=true)
        dist_points = points_start[:, i_sample]
        dist_to_centre = compute_dist_to_centre(point_centre, dist_points, R)
        p = density(dist_to_centre, bandwidth=0.05)

        data_save = DataFrame(Dict("dist_centre" => dist_to_centre))
        path_save = joinpath(dir_save, "Dist_data_weight$power_weight.csv")
        CSV.write(path_save, data_save)
    end
end

if ARGS[1] == "Illustration_particle_diffusion_recruited" || ARGS[1] == "all_illustrations"
    #=
        Simulate particle model with diffusion
    =#
    # Parameter needed by simulator
    dir_save = joinpath(@__DIR__, "..", "..", "Results", "Particle_simulator_diff", "Illustration_recruited")
    if !isdir(dir_save)
        mkpath(dir_save)
    end
    Δt = 0.01
    frac_cover=0.01
    n_particles_deliver=100
    k_off = 0.1
    power_weight = 3.0
    t_end = 400.0
    λhit=0.4
    n_save=25
    R=2.5
    α=0.2
    γ=1.0
    have_exo=true
    recruit_with_exo = true
    recruit_particles_reaction = true

    k_off = 0.5
    k_recruit = 0.0
    data_save = Vector{Vector{Matrix{Float64}}}(undef, 2)

    # Fast Dm
    Dm = 0.0045
    Random.seed!(123)
    res = simulate_particles(Δt, frac_cover, n_particles_deliver, k_off, power_weight,
                            t_end, λhit, n_save, Dm, R, α, γ, k_recruit, have_exo, recruit_with_exo, recruit_particles_reaction, false)

    i_list = [1, 2, 4, 5, 10]
    for i in i_list
        _t = (t_end - 0.0) / n_save * (i - 1)
        dist_to_centre = compute_dist_to_centre(point_centre, res[i], R)
        p1 = surface(sphere(2.5, [0.0, 0.0, 0.0]),
                    camera=(45, 10),
                    xlim=[-R, R].*1.2,
                    ylim=[-R, R].*1.2,
                    zlim=[-R, R].*1.2,
                    fill="grey",
                    xlab="x", ylab="y", zlab="z",
                    title = "t = $_t [s]",
                    legend=false)
        p1 = scatter3d!(res[i][1, :], res[i][2, :], res[i][3, :],
                        alpha=0.7,
                        camera=(45, 10))
        path_save = joinpath(dir_save, "t$_t.svg")
        savefig(p1, path_save)

        p2 = density(dist_to_centre, title = "t = $_t", label=false)
        path_save = joinpath(dir_save, "Dens_t$_t.svg")
        savefig(p2, path_save)
    end
end
