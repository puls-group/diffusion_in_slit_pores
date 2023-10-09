#!/usr/bin/env python3
import scipy.integrate as sp_int
import numpy as np
from os.path import exists, isfile
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from matplotlib import cm
from scipy.interpolate import interp1d

cache_limit = 1.1
loaded_cache_data = None
rel_displ = None
rel_diff = None
rel_lt = None

local_cache_path = "./cache/local_lv_correction.cache"


"""
TODO: replace with LV calculation code

def simulate_fpe_for(
    rel_displacement: float,
    rel_D_displacement: float,
    dx=0.01,
    dt=0.1,
    return_detailed_stats=False,
):
    # We set the slab width to 1 to simplify the model, everything scales relatively, so that's fine
    slab_width = 1

    # Basic dimension of the domain in z direction and d direction accounting for dx resolution
    dim_z = int(np.ceil(slab_width / dx))
    dim_d = int(np.ceil(rel_displacement / dx))

    # Total dimension of the domain in d direction accounting for both signs of displacement
    dim_s = 2 * dim_d + 1
    # Total dimension in z direction accounting for tilted shape
    total_z = dim_z + dim_s + 2
    total_entry = dim_s * total_z

    # Set the diffusion parameters for the simulation
    D_translational = 1
    D_displacement = rel_D_displacement

    lt_translational_estimate = 1.0 / 12.0 * slab_width ** 2 / D_translational

    print(D_translational, D_displacement, slab_width, rel_displacement)

    domain = np.zeros((dim_s, total_z))
    mask = domain * 0.0

    D_translational /= dx * dx
    D_displacement /= dx * dx

    dt = min(dt, lt_translational_estimate / 100.0)

    total_norm = 0.0
    # initialize domain
    for s in range(dim_s):
        start_off = s + 1
        domain[s, start_off: (start_off + dim_z)] = 1.0
        mask[s, start_off: (start_off + dim_z)] = 1.0
        total_norm += dim_z

    def fpe_diff(t, cy):
        # global D_translational, D_displacement, dim_s, total_z, total_entry
        y = cy.reshape((dim_s, total_z))
        diff = y * 0.0

        s_up = y[2:]
        s_mid = y[1:-1]
        s_down = y[:-2]
        diff_s_mid = s_up - 2 * s_mid + s_down
        diff[1:-1] += D_displacement * diff_s_mid

        # FIXME: No factor 2 because there is only flow in one of two possible directions?
        # diff[0] += 2.0 * D_displacement * (y[1] - y[0])
        diff[0] += D_displacement * (y[1] - y[0])
        # diff[-1] += 2.0 * D_mol * (y[-2] - y[-1])
        diff[-1] += D_displacement * (y[-2] - y[-1])

        z_up = y[:, 2:]
        z_mid = y[:, 1:-1]
        z_down = y[:, :-2]
        diff[:, 1:-1] += D_translational * (z_up - 2 * z_mid + z_down)

        diff = diff * mask
        return diff.reshape((total_entry,))

    # print(total_norm, np.sum(domain))
    # sys.exit(1)

    t0 = 0.0
    y0 = domain.reshape((total_entry,))

    t_bound = max(5.0, lt_translational_estimate * 8)
    max_step = dt

    print(lt_translational_estimate)
    # sys.exit()

    last_report = 0.0
    report_step = t_bound / 100.0

    sim_lt = [t0]
    survival_probability = [1.0]
    total_lifetime_integral = 0.0

    FPE = sp_int.RK45(fpe_diff, t0, y0, t_bound=t_bound, max_step=max_step)
    # pprint(FPE.__dict__)
    print(FPE.status)

    while True:
        status_update = FPE.step()
        # pprint(FPE.__dict__)
        curr_p = np.sum(FPE.y) / total_norm

        total_lifetime_integral += (
            (survival_probability[-1] + curr_p) * (FPE.t - sim_lt[-1]) / 2.0
        )

        sim_lt.append(FPE.t)
        survival_probability.append(curr_p)

        # print("{0:.2f}/{1:.2f}".format(FPE.t, t_bound))

        if last_report + report_step <= FPE.t:
            domain = FPE.y.reshape((dim_s, total_z))

            last_report += report_step

        if FPE.status != "running" or curr_p < 1e-6:
            print(FPE.status, status_update)
            break

    est_total_time = total_lifetime_integral
    # Calculate exponential approximation integral for infinite tail
    num_entries = len(sim_lt)

    approx_entries = int(num_entries * 2.0 / 3.0)
    coeff = np.polyfit(
        sim_lt[approx_entries:],
        np.log(survival_probability[approx_entries:]),
        deg=1,
    )
    f_lin = np.poly1d(coeff)

    est_total_time += np.exp(f_lin(sim_lt[-1])) / coeff[1]
    rel_total_time = est_total_time / lt_translational_estimate

    print(est_total_time, "/", lt_translational_estimate, "=", rel_total_time)

    if return_detailed_stats:
        return (
            rel_total_time,
            FPE.t / lt_translational_estimate,
            (sim_lt, survival_probability),
        )
    else:
        return (rel_total_time, FPE.t / lt_translational_estimate)

"""

def correction(
    D_trans: float,
    slab_width: float,
    D_displacement: float,
    displacement: float,
    skip_cache: bool = False,
) -> float:
    global local_cache_path
    rel_displacement = displacement / slab_width
    rel_diffusion = D_displacement / D_trans

    global cache_limit, loaded_cache_data, rel_displ, rel_diff, rel_lt

    if exists(local_cache_path) and isfile(local_cache_path):
        if not skip_cache:
            if loaded_cache_data is None:
                cached_data = np.loadtxt(local_cache_path, ndmin=2)
                if len(cached_data) > 0:
                    rel_displ = cached_data[:, 0]
                    rel_diff = cached_data[:, 1]
                    rel_lt = cached_data[:, 2]

                    rel_disp_filter = (rel_displ <= cache_limit)

                    filter_data = cached_data[rel_disp_filter]

                    rel_displ = filter_data[:, 0]
                    rel_diff = filter_data[:, 1]
                    rel_lt = filter_data[:, 2]

                    cached_displ_set = list(dict.fromkeys(rel_displ))

                    cache_inter_func = []
                    cache_inter_displ = []

                    # Build interpolation functions:

                    for i in range(len(cached_displ_set)):
                        picked_displ = cached_displ_set[i]

                        picked_filter = (rel_displ > (picked_displ -
                                                      1e-2)) & (rel_displ < (picked_displ + 1e-2))

                        picked_data = filter_data[picked_filter]

                        # We need a minimum of values to operate
                        if len(picked_data) < 3:
                            continue

                        picked_data = picked_data[picked_data[:, 1].argsort(
                        )]

                        picked_diff = picked_data[:, 1]
                        picked_lt = picked_data[:, 2]
                        picked_lt = picked_lt/max(1.0, np.max(picked_lt))

                        limit = (
                            1.-2.*picked_displ)**3 if picked_displ < 0.5 else 0.0

                        # We do relative interpolation between 1 and (1-2r)^3 limits
                        picked_weight = (picked_lt-limit)/(1.-limit)

                        interpolator = interp1d(
                            picked_diff, picked_weight, fill_value=(1.0, 0.0), bounds_error=False)

                        upper_value_index = np.argmax(picked_diff)
                        upper_diff = picked_diff[upper_value_index]
                        upper_weight = (
                            picked_lt[upper_value_index]-limit)/(1.-limit)

                        def f(x):
                            if x < upper_diff:
                                return interpolator(x)
                            else:
                                return (upper_weight * upper_diff / x)

                        cache_inter_func.append(f)
                        cache_inter_displ.append(cached_displ_set[i])

                def apply_cache_data(rel_d, rel_D):
                    ref_limit = (1.-2.*rel_d)**3 if rel_d < 0.5 else 0.0

                    # Deal with low values
                    if rel_d < cached_displ_set[0]:
                        interp_weight = (
                            rel_d * cache_inter_func[0](rel_D) + (cached_displ_set[0]-rel_d) *1.)/cached_displ_set[0]
                        return (interp_weight * (1.-ref_limit))+ref_limit

                    # Deal with upper convergence
                    if rel_d >= cached_displ_set[-1]:
                        interp_weight = cache_inter_func[-1](rel_D)
                        return interp_weight

                    # Interpolate in between
                    for i in range(len(cached_displ_set)-1):
                        if rel_d >= cached_displ_set[i] and rel_d < cached_displ_set[i+1]:
                            diff_displ = cached_displ_set[i+1] - \
                                cached_displ_set[i]
                            # We interpolate linearly for the weight in between 1 and (1-2r)**3:
                            interp_weight = ((cached_displ_set[i+1]-rel_d) *
                                             cache_inter_func[i](
                                rel_D) + (rel_d-cached_displ_set[i]) *
                                cache_inter_func[i+1](
                                    rel_D))/diff_displ
                            # Then calculate the actual approximation from the interpolation
                            return (interp_weight * (1.-ref_limit))+ref_limit
                loaded_cache_data = apply_cache_data

            # if rel_displacement <= cache_limit:
            return loaded_cache_data(rel_displacement, rel_diffusion)

    else:
        with open(local_cache_path, "w") as cache_out:
            cache_out.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                    "#relative displacement",
                    "#relative deform diffusion",
                    "#relative lifetime compared to expected lifetime",
                    "#simulated time relative to expected lifetime",
                    "#time spent during simulation [s]",
                )
            )

    start = timer()
    total_rel_time, total_sim_limit = simulate_fpe_for(
        rel_displacement, rel_diffusion)
    end = timer()

    sim_time = end - start
    print("Spent", sim_time, "during simulation")

    with open(local_cache_path, "a") as cache_out:
        cache_out.write(
            "{0:.8e}\t{1:.8e}\t{2:.8e}\t{3:.8e}\t{4:.8e}\n".format(
                rel_displacement,
                rel_diffusion,
                total_rel_time,
                total_sim_limit,
                sim_time,
            )
        )

    return total_rel_time

