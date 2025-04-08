import uproot
import matplotlib.pyplot as plt
import numpy
import math
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.lines import Line2D

# crossection = 3714 * 1000
# Ngen = 99500
# Ldat = 1.7
# w2 = Ldat / (Ngen / crossection)
# file_uproot = uproot.open("output-egamma.root")

# divide_data_pt = file_uproot["d_rec_eff_data_pt"]
# divide_data_eta = file_uproot["d_rec_eff_data_eta"]
# divide_mc_pt = file_uproot["d_rec_eff_mc_pt"]
# divide_mc_eta = file_uproot["d_rec_eff_mc_eta"]
# d_ratio_pt = file_uproot["d_ratio_pt"]
# d_ratio_eta = file_uproot["d_ratio_eta"]

# wils_rec_mc_pt = file_uproot["d_rec_mc_pt"]
# wils_rec_mc_eta = file_uproot["d_rec_mc_eta"]
# wils_rec_data_pt = file_uproot["d_rec_data_pt"]
# wils_rec_data_eta = file_uproot["d_rec_data_eta"]

# wils_probe_mc_pt = file_uproot["d_probe_mc_pt"]
# wils_probe_mc_eta = file_uproot["d_probe_mc_eta"]
# wils_probe_data_pt = file_uproot["d_probe_data_pt"]
# wils_probe_data_eta = file_uproot["d_probe_data_eta"]

# op_vl = file_uproot["d_rec_eff_data_pt_vl"]
# op_lo = file_uproot["d_rec_eff_data_pt_lo"]
# op_me = file_uproot["d_rec_eff_data_pt_me"]
# op_ti = file_uproot["d_rec_eff_data_pt_ti"]


def ratio_1d(hist_name, xtit, fname):
    plt.rcParams.update({'font.size': 7})
    files = ["output-basev2.root", "output-acov2.root", "output-pixv2.root", "output-trigdownv2.root", "output-trigupv2.root"]
    
    # Load the base histogram and retrieve data
    file = uproot.open("output-basev2.root")
    hist = file[hist_name]
    contents, edges = hist.to_numpy()
    variances = hist.variances()
    uncertainties = np.sqrt(variances)
    # contents = contents[2:]
    # edges = edges[2:]
    # uncertainties = uncertainties[2:]
    
    # Initialize uncertainty arrays for upper and lower bounds
    uncertainties_bot = np.zeros_like(uncertainties)
    uncertainties_top = np.zeros_like(uncertainties)
    err_vec_top, err_vec_bot = syst_err_1d(files, hist_name)
    # err_vec_bot = err_vec_bot[2:]
    # err_vec_top = err_vec_top[2:]
    fig, ax = plt.subplots(dpi=300)
    
    # Calculate total uncertainties by combining systematic and statistical uncertainties
    for i in range(len(contents)):
        uncertainties_top[i] = math.sqrt(err_vec_top[i] ** 2 + uncertainties[i] ** 2)
        uncertainties_bot[i] = math.sqrt(err_vec_bot[i] ** 2 + uncertainties[i] ** 2)
    
    x_centers = (edges[:-1] + edges[1:]) / 2
    x_unc = (edges[1:] - edges[:-1]) / 2

    # Plot data points with error bars
    ax.errorbar(
        x_centers, contents, xerr=x_unc,
        yerr=[uncertainties_bot, uncertainties_top], 
        fmt='o', color='black', ecolor='black', capsize=1,
    )

    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1)
    print(edges[0])
    #Annotate each bin with content and uncertainties
    for i in range(len(contents)):
        offset = 0.05
        # if i==1 or i==8:
        #     offset = -0.3
        x_center = (edges[i] + edges[i + 1]) / 2
        content = contents[i]
        uncertainty_top = uncertainties_top[i]
        uncertainty_bot = uncertainties_bot[i]
        # Display the content value at the point
        ax.text(x_center, content - uncertainty_bot - 0.05, f"{content:.3f}", ha="center", va="center", color="red", size=5., weight='bold')
        # Display the positive uncertainty above the point
        # ax.text(x_center, content + uncertainty_top+0.005, f"+{uncertainty_top:.2f}", 
        #         ha="center", va="bottom", color="blue", size=4.)
        
        # # Display the negative uncertainty below the point
        # ax.text(x_center, content - uncertainty_bot-0.005, f"-{uncertainty_bot:.2f}", 
        #         ha="center", va="top", color="blue", size=4.)


    # Configure plot scale and labels
    ax.set_ylabel("Scale Factor")
    ax.set_xlabel(xtit)
    # ax.set_xlim(0., 50.)
    # ax.set_ylim(0.8, 1.3)
    print(fname)
    plt.savefig(fname+".png")
    plt.show()

def syst_err_1d(files, hist_name):
    files_len = len(files)
    mats = []
    
    # Load each file and get 1D histogram data
    for file in files:
        hist = uproot.open(file)[hist_name]
        contents, _ = hist.to_numpy()
        print(contents.shape)
        mats.append(contents)
    
    mats = np.array(mats)
    size = mats.shape[1]
    
    # Calculate systematic uncertainties (upper and lower) for each bin
    err_vec_bot = np.zeros(size)
    err_vec_top = np.zeros(size)
    
    for i in range(size):
        for k in range(1, files_len):
            diff = mats[0][i] - mats[k][i]
            if diff < 0:
                err_vec_top[i] += diff ** 2
            else:
                err_vec_bot[i] += diff ** 2
    
    err_vec_bot = np.sqrt(err_vec_bot)
    err_vec_top = np.sqrt(err_vec_top)
    return err_vec_top, err_vec_bot


#ratio_1d("d_ratio_pt", "$p_T^{probe} [GeV]$", "scale_factor_pt")
ratio_1d("d_ratio_eta", "$\eta^{probe}$", "scale_factor_eta")


def ratio_2d():
    plt.rcParams.update({'font.size': 5})
    files = ["output-basev2.root", "output-acov2.root", "output-pixv2.root", "output-trigdownv2.root", "output-trigupv2.root"]
    file = uproot.open("output-basev2.root")
    hist = file["d_ratio_2d"]
    contents, x_edges, y_edges = hist.to_numpy()
    variances = hist.variances()
    uncertainties = np.sqrt(variances)
    size_unc = uncertainties.shape
    uncertainties_bot = numpy.zeros(size_unc)
    uncertainties_top = numpy.zeros(size_unc)
    err_mat_top, err_mat_bot = syst_err(files)
    fig, ax = plt.subplots(dpi=200)
    for i in range(len(contents)):
        for j in range(len(contents[0])):
            uncertainties_top[i, j] = math.sqrt(err_mat_top[i, j] ** 2 + uncertainties[i, j] ** 2)
            uncertainties_bot[i, j] = math.sqrt(err_mat_bot[i, j] ** 2 + uncertainties[i, j] ** 2)
    contents[0,:] = np.NAN
    contents[:,1] = np.NAN
    contents[:,8] = np.NaN
    mesh = ax.pcolormesh(x_edges, y_edges, contents.T, cmap="Greys", shading="auto")
    print(err_mat_top)
    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label("Electron reconstruction efficiency")

    for i in range(len(x_edges) - 1):
        for j in range(len(y_edges) - 1):
            if not((i==0) or (j==8 or j==1)):
                x_center = numpy.exp((numpy.log(x_edges[i]) + numpy.log(x_edges[i + 1])) / 2)
                y_center = (y_edges[j] + y_edges[j + 1]) / 2
                content = contents[i, j]
                uncertainty_top = uncertainties_top[i, j]
                uncertainty_bot = uncertainties_bot[i ,j]
                ax.text(x_center, y_center, f"{content:.3f}\n+{uncertainty_top:.3f}\n-{uncertainty_bot:.3f}", 
                        ha="center", va="center", color="Red", size=5., weight='bold')
    ax.set_xscale('log')
    tick_positions = [0.5, 2, 3, 5, 7, 10, 50]
    tick_labels = ['0.5', '2', '3', '5', '7', '10', '50']
    ax.xaxis.set_major_locator(FixedLocator(tick_positions))
    ax.xaxis.set_major_formatter(FixedFormatter(tick_labels))
    tick_positions = [-2.47, -1.52, -1.37 , -1.,  -0.5, 0., 0.5, 1., 1.37, 1.52, 2.47]
    tick_labels = ["-2.47", "-1.52", "-1.37 ", "-1.", " -0.5", "0.", "0.5", "1.", "1.37", "1.52", "2.4"]
    ax.yaxis.set_major_locator(FixedLocator(tick_positions))
    ax.yaxis.set_major_formatter(FixedFormatter(tick_labels))
    
    # Add gridlines to separate cells
    ax.grid(which='major', color='black', linestyle='--', linewidth=0.5, alpha=0.5)

    # Ensure the grid aligns with the edges of the mesh
    ax.set_xlim(2., x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])
    ax.set_xlabel("$p_{T}^{probe}$ [GeV]")
    ax.set_ylabel("$\eta^{probe}$")
    plt.savefig("scale_factor_2d.png")
    plt.show()

def syst_err(files):
    files_len = len(files)
    files_uproot = [None] * files_len
    hists = [None] * files_len
    mat = [None] * files_len
    for i in range(files_len):
        files_uproot[i] = uproot.open(files[i])
        hists[i] = files_uproot[i]["d_ratio_2d"]
        mat[i],_ ,_ = hists[i].to_numpy()
    mat = numpy.array(mat)
    size = mat[0].shape
    err_mat_bot = numpy.zeros(size)
    err_mat_top = numpy.zeros(size)
    for i in range(len(mat[0])):
        for j in range(len(mat[0][0])):
            for k in range(1, len(mat)):
                mat[k][i][j] = mat[0][i][j] - mat[k][i][j]
                if(mat[k][i][j] < 0):
                    err_mat_top[i][j] += mat[k][i][j] ** 2
                else:
                    err_mat_bot[i][j] += mat[k][i][j] ** 2
    err_mat_bot = numpy.sqrt(err_mat_bot)
    err_mat_top = numpy.sqrt(err_mat_top)
    return err_mat_top, err_mat_bot

#ratio_2d()

def calculate_efficiency_wilson(num, den):
    # Number of bins
    n_bins = len(num)

    # Initialize arrays for efficiency and errors
    eff = np.zeros(n_bins)
    eff_err = np.zeros(n_bins)

    # Set z-value for 1-sigma confidence level (68%)
    z = 1.0  # 1-sigma (68%) corresponds to z = 1.0

    for i in range(n_bins):
        n = num[i]  # Numerator
        d = den[i]  # Denominator

        if d == 0:  # Avoid division by zero
            eff[i] = 0
            eff_err[i] = 0
        else:
            phat = n / d  # Efficiency
            eff[i] = phat

            # Wilson score interval calculation with 1-sigma confidence (z = 1.0)
            denominator = 1 + z*z/d
            center = (phat + z*z/(2*d)) / denominator
            margin = z * np.sqrt((phat*(1-phat) + z*z/(4*d)) / d) / denominator

            # Calculate bounds and symmetric error approximation
            lower_bound = center - margin
            upper_bound = center + margin

            # Symmetric error approximation (max difference from center)
            eff_err[i] = max(upper_bound - phat, phat - lower_bound)

    return eff, eff_err


def plot_eff(data, mc, xlabel, tit, files, hist_name_mc, hist_name_data):
    # Extract bin edges, values, and statistical errors
    edges = data.axis().edges()
    data_val = data.values()
    data_stat_err = np.sqrt(data.variances())
    mc_val = mc.values()
    mc_stat_err = np.sqrt(mc.variances())
    
    # Calculate systematic errors
    data_syst_err_top, data_syst_err_bot = syst_err_1d(files, hist_name_data)
    mc_syst_err_top, mc_syst_err_bot = syst_err_1d(files, hist_name_mc)
    
    # Combine statistical and systematic uncertainties
    data_err_top = np.sqrt(data_stat_err**2 + data_syst_err_top**2)
    data_err_bot = np.sqrt(data_stat_err**2 + data_syst_err_bot**2)
    mc_err_top = np.sqrt(mc_stat_err**2 + mc_syst_err_top**2)
    mc_err_bot = np.sqrt(mc_stat_err**2 + mc_syst_err_bot**2)
    
    # Calculate bin centers and bin width (for xerror)
    x = (edges[:-1] + edges[1:]) / 2.0  # Update x-bin centers after removing 2 bins
    xerror = (edges[1:] - edges[:-1]) / 2.0  # Update x-error (half-width of the bins)

    # Create the figure and subplots, sharing the X-axis
    fig, axs = plt.subplots(dpi=200)

    # Plot Data efficiency with combined errors
    axs.errorbar(
        x, data_val, xerr=xerror, yerr=[data_err_bot, data_err_top], 
        fmt='o', color='black', ecolor='black', capsize=2, label="Data"
    )

    # Plot MC efficiency with combined errors
    axs.errorbar(
        x, mc_val, xerr=xerror, yerr=[mc_err_bot, mc_err_top], 
        fmt='o', color='blue', ecolor='blue', capsize=2, label="MC"
    )

    # Draw a line at y=1 for reference
    axs.axhline(y=1, color='gray', linestyle='--', linewidth=1)

    # Set labels, title, and legend
    axs.set_xlabel(xlabel)
    axs.set_ylabel("Electron reconstruction efficiency")
    axs.legend(loc='best')

    # Set y-axis limits for the plot
    axs.set_ylim(0, 1)
    # Adjust layout to avoid overlap
    plt.tight_layout()
    plt.savefig(tit)
    plt.show()


# file1 = uproot.open("output-basev2.root")
# files = ["output-basev2.root", "output-acov2.root", "output-pixv2.root", "output-trigdownv2.root", "output-trigupv2.root"]
# plot_eff(file1["d_rec_eff_data_pt"], file1["d_rec_eff_mc_pt"], "$p_T^{probe} [GeV]$", "efs_pt", files, "d_rec_eff_mc_pt", "d_rec_eff_data_pt")
# plot_eff(file1["d_rec_op_data_eta"], file1["d_rec_eff_mc_eta"], "$\eta^{probe}$", "efs_eta", files, "d_rec_eff_mc_eta", "d_rec_op_data_eta")



def plot_wils(data_n, data_d, mc_n, mc_d, xlabel, tit):

    data_val, data_err = calculate_efficiency_wilson(data_n.values(), data_d.values())
    mc_val, mc_err = calculate_efficiency_wilson(mc_n.values(), mc_d.values())
    ratio_val, ratio_err = calculate_ratio(data_val, data_err, mc_val, mc_err)

    edges = data_n.axis().edges()

    # Create the figure and subplots, sharing the X-axis
    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(6, 8), sharex=True)

    # Calculate bin centers and bin width (for xerror)
    x = (edges[:-1] + edges[1:]) / 2.0
    xerror = (edges[1:] - edges[:-1]) / 2.0

    # Upper plot: Data and MC efficiency
    axs[0].errorbar(x, data_val, xerr=xerror, yerr=data_err, fmt='o', color='black', ecolor='black', label="Data")
    axs[0].errorbar(x, mc_val, xerr=xerror, yerr=mc_err, fmt='o', color='blue', ecolor='blue', label="MC")
    axs[0].legend(loc='best')
    axs[0].set_ylabel("Electron reconstruction efficiency")
    #axs[0].text(0.95, 0.05, "ATLAS", verticalalignment='bottom', horizontalalignment='right',
    #            transform=axs[0].transAxes, fontsize=10, color='black')


    # Set y-axis limits for the upper plot
    axs[0].set_ylim(0, 1)

    # Lower plot: Data/MC ratio
    axs[1].errorbar(x, ratio_val, xerr=xerror, yerr=ratio_err, fmt='o', color='black', ecolor='black')
    axs[1].set_ylabel("Data/MC")

    # Set y-axis limits for the lower plot
    axs[1].set_ylim(0.7, 1.3)

    # Add a horizontal red line at y = 1.0 in the lower plot
    axs[1].axhline(1.0, color='red', linestyle='--', linewidth=1)

    # Common X-axis label under the lower plot
    axs[1].set_xlabel(xlabel)

    # Adjust layout to avoid overlap
    plt.tight_layout()
    plt.savefig(tit)
    plt.show()


def calculate_ratio(num_eff, num_err, den_eff, den_err):
    # Number of bins
    n_bins = len(num_eff)
    
    # Initialize arrays for ratio and ratio error
    ratio = np.zeros(n_bins)
    ratio_err = np.zeros(n_bins)

    for i in range(n_bins):
        if den_eff[i] == 0:
            # If denominator efficiency is zero, set ratio and error to 0
            ratio[i] = 0
            ratio_err[i] = 0
        else:
            # Calculate the ratio
            ratio[i] = num_eff[i] / den_eff[i]

            # Propagate the error using the formula
            ratio_err[i] = ratio[i] * np.sqrt(
                (num_err[i] / num_eff[i])**2 + (den_err[i] / den_eff[i])**2
            )

    return ratio, ratio_err


# plot_eff(data=divide_data_pt, mc=divide_mc_pt, ratio=d_ratio_pt, xlabel="electron $p^{probe}_{T}$ [GeV]", tit="divide_pt.png")
# plot_eff(data=divide_data_eta, mc=divide_mc_eta, ratio=d_ratio_eta, xlabel="electron $\eta^{probe}$", tit="divide_eta.png")
# plot_wils(data_n=wils_rec_data_pt, data_d=wils_probe_data_pt, mc_n=wils_rec_mc_pt, mc_d=wils_probe_mc_pt, xlabel="electron $p^{probe}_{T}$ [GeV]", tit="wilson_pt.png")
# plot_wils(data_n=wils_rec_data_eta, data_d=wils_probe_data_eta, mc_n=wils_rec_mc_eta, mc_d=wils_probe_mc_eta, xlabel="electron $\eta^{probe}$", tit="wilson_eta.png")



# mc_val = wils_probe_mc_pt.values()
# data_val = wils_probe_data_pt.values()
# edges = wils_probe_data_pt.axis().edges()
# plt.figure()
# mc_val = w2 * mc_val
# plt.hist(edges[:-1], bins=edges, weights=mc_val, color='blue', histtype='step', label='MC')
# plt.hist(edges[:-1], bins=edges, weights=data_val, color='black', histtype='step', label='Data')
# plt.xlabel("electron $p_{T}^{probe}$ [GeV]")
# plt.ylabel("No electrons")
# plt.legend(loc='best')
# plt.yscale('log')
# plt.savefig("NoProbes.png")
# plt.show()


# vl = op_vl.values()
# lo = op_lo.values()
# me = op_me.values()
# ti = op_ti.values()
# edges = op_vl.axis().edges()
# vl_err = np.sqrt(op_vl.variances())
# lo_err = np.sqrt(op_lo.variances())
# me_err = np.sqrt(op_me.variances())
# ti_err = np.sqrt(op_ti.variances())
# x = (edges[:-1] + edges[1:]) / 2.0
# xerror = (edges[1:] - edges[:-1]) / 2.0
# plt.errorbar(x, vl, xerr=xerror, yerr=vl_err, label="Very Loose", fmt='o')
# plt.errorbar(x, lo, xerr=xerror, yerr=lo_err, label="Loose", fmt='o')
# plt.errorbar(x, me, xerr=xerror, yerr=me_err, label="Medium", fmt='o')
# plt.errorbar(x, ti, xerr=xerror, yerr=ti_err, label="Tight", fmt='o')
# plt.legend(loc='best')
# plt.ylabel("Electron identification efficiency")
# plt.ylim(0, 1)
# plt.xlabel("electron $p^{probe}_{T}$ [GeV]")
# plt.savefig("op.png")
# plt.show()

# aco_truth = file_uproot["d_tp_mc_aco"]
# aco_data = file_uproot["d_tp_data_aco"]
# aco_cutvl = file_uproot["d_cut_data_aco"]
# tclu_mult = file_uproot["d_tclu_mult"]
# tclu_mult_nvl = file_uproot["d_tclu_mult_nvl"]

# aco_truth_val = aco_truth.values()
# aco_data_val = aco_data.values()
# aco_cutvl_val = aco_cutvl.values()
# edges = aco_truth.axis().edges()

# aco_truth_val /= np.sum(aco_truth_val)
# aco_data_val /= np.sum(aco_data_val)
# aco_cutvl_val /= np.sum(aco_cutvl_val)

# plt.figure(dpi=200)
# plt.hist(edges[:-1], bins=edges, weights=aco_data_val, color='black', histtype='step', label='Data')
# plt.hist(edges[:-1], bins=edges, weights=aco_truth_val, color='red', histtype='step', label='MC')
# plt.hist(edges[:-1], bins=edges, weights=aco_cutvl_val, color='blue', histtype='step', label='Data, cutted by VL')
# plt.legend(loc='best')
# plt.ylabel("Events")
# plt.yscale('log')
# plt.xlabel("Aco tag-probe")
# plt.savefig("aco_comp.png")
# plt.show()

# tclu_mult_val = tclu_mult.values()
# tclu_mult_nvl_val = tclu_mult_nvl.values()
# edges = tclu_mult.axis().edges()
# plt.figure(dpi=200)
# tclu_mult_val /= np.sum(tclu_mult_val)
# tclu_mult_nvl_val /= np.sum(tclu_mult_nvl_val)
# plt.hist(edges[:-1], bins=edges, weights=tclu_mult_val, histtype='step', label="Tag-probe")
# plt.hist(edges[:-1], bins=edges, weights=tclu_mult_nvl_val, histtype='step', label="Not passed vl")
# plt.legend(loc='best')
# plt.ylabel("Events")
# plt.yscale('log')
# plt.xlabel("No Topoclusters")
# plt.savefig("tclu_mult.png")
# plt.show()



# file = uproot.open("output-basev2.root")
# hist = file["d_ratio_2d"]  # Replace "hist_name" with your histogram's name
# contents, x_edges, y_edges = hist.to_numpy()

# # Create a DataFrame
# data = []

# # Iterate over all bins
# for i in range(len(x_edges) - 1):  # Loop over x bins
#     for j in range(len(y_edges) - 1):  # Loop over y bins
#         data.append({
#             "X Bin Start": x_edges[i],
#             "X Bin End": x_edges[i + 1],
#             "Y Bin Start": y_edges[j],
#             "Y Bin End": y_edges[j + 1],
#             "Content": contents[i, j]
#         })

# # Convert to DataFrame
# df = pd.DataFrame(data)

# # Save to a CSV file
# df.to_csv("2D_scale_factor.csv", index=False)

# file = uproot.open("output-basev2.root")
# hist = file["d_ratio_pt"]  # Replace "hist_name" with your histogram's name
# contents, edges = hist.to_numpy()

# # Create a DataFrame
# df = pd.DataFrame({
#     "Bin Start": edges[:-1],  # Left edges of bins
#     "Bin End": edges[1:],    # Right edges of bins
#     "Content": contents      # Bin contents
# })

# # Save to a CSV file
# df.to_csv("scale_factor_pt.csv", index=False)

# hist = file["d_ratio_eta"]  # Replace "hist_name" with your histogram's name
# contents, edges = hist.to_numpy()

# # Create a DataFrame
# df = pd.DataFrame({
#     "Bin Start": edges[:-1],  # Left edges of bins
#     "Bin End": edges[1:],    # Right edges of bins
#     "Content": contents      # Bin contents
# })

# # Save to a CSV file
# df.to_csv("scale_factor_eta.csv", index=False)


def ratio_1d(hist_name, xtit, fname):
    plt.rcParams.update({'font.size': 5})
    values_arr = [None] * 7
    files = numpy.array(["output-basev2.root", "output-acov2.root", "output-pixv2.root", "output-tigv2.root", "output-trigdownv2.root", "output-trigupv2.root", "output-zdcv2.root"])
    plt.figure(dpi=200)
    for i in range(files.shape[0]):
        filename = files[i]
        file_uproot = uproot.open(filename)
        hist = file_uproot[hist_name]
        contents, edges = hist.to_numpy()
        x_centers = (edges[:-1] + edges[1:]) / 2
        values_arr[i] = contents
        if i:
            values_arr[i] = numpy.abs(values_arr[i]-values_arr[0]) / values_arr[0]
    file_uproot = uproot.open("output-basev2.root")
    hist = file_uproot[hist_name]
    values_arr[0] = hist.variances()
    values_arr = numpy.array(values_arr)
    # values_arr = np.delete(values_arr, 1, axis=1)
    # values_arr = np.delete(values_arr, 8, axis=1)
    # x_centers = np.delete(x_centers,1)
    # x_centers = np.delete(x_centers,7)
    # print(x_centers)
    # x_centers = x_centers[2:]
    # values_arr = values_arr[:,2:]
    for i in range(values_arr.shape[0]):
        plt.plot(x_centers, values_arr[i])
    plt.legend(["Stat", "Aco", "Pixtrack", "LHTight", "Trigger_down", "Trigger_up", "ZDC"], ncol=2, loc="best")
    plt.xlabel(xtit)
    plt.ylabel("Relative Uncertainty")
    plt.savefig(fname+".png")
    plt.show()


    # values_arr = numpy.array(values_arr)
    # df = pd.DataFrame({
    # "Bin Start": edges[:-1],  # Left edges of bins
    # "Bin End": edges[1:],    # Right edges of bins
    # "Base": values_arr[1,:],   
    # "Aco": values_arr[2,:],
    # "Pix": values_arr[3,:],
    # "Tig": values_arr[4,:],
    # "TrigD": values_arr[5,:],
    # "TrigU": values_arr[6,:]
    # })
    # df.to_csv(fname+".csv", index=False)
    

# ratio_1d("d_ratio_pt", "$p_T^{probe} [GeV]$", "pt_sf_com")
# ratio_1d("d_ratio_eta", "$\eta^{probe} $", "eta_sf_com")

def control_plots(histname_mc, histname_data, filename, xaxis_label, xlim):
    files = ["output-basev2.root", "output-tigv2.root"]
    contents_mc = []
    edges_mc = []
    for i in range(2):
        files_uproot = uproot.open(files[i])
        hist = files_uproot[histname_mc]
        content, edge = hist.to_numpy()
        contents_mc.append(content)
        edges_mc.append(edge)
    contents_mc = numpy.array(contents_mc)
    edges_mc = numpy.array(edges_mc)
    contents_data = []
    edges_data = []
    for i in range(2):
        files_uproot = uproot.open(files[i])
        hist = files_uproot[histname_data]
        content, edge = hist.to_numpy()
        contents_data.append(content)
        edges_data.append(edge)
    for i in range(2):
        contents_mc[i] /= numpy.sum(contents_mc[i])
        contents_data[i] /= numpy.sum(contents_data[i])
    contents_data = numpy.array(contents_data)
    edges_data = numpy.array(edges_data)
    x_centers = (edges_mc[0][:-1] + edges_mc[0][1:]) / 2

    # Main plots
    fig, ax = plt.subplots(2, 1, figsize=(8, 10), dpi=150)
    ax[0].errorbar(x_centers, contents_data[0], fmt='o', color='black', label="Data", markersize=4)
    ax[0].hist(edges_mc[0][:-1], bins=edges_mc[0], weights=contents_mc[0], color='red', histtype='step', label="MC")
    ax[0].set_ylabel("Events Medium")
    ax[0].set_yscale('log')
    ax[0].set_xlabel(xaxis_label)
    ax[0].legend()
    ax[0].set_xlim(xlim)
    ax[1].errorbar(x_centers, contents_data[1], fmt='o', color='black', label="Data", markersize=4)
    ax[1].hist(edges_mc[0][:-1], bins=edges_mc[0], weights=contents_mc[1], color='red', histtype='step', label="MC")
    ax[1].set_ylabel("Events Tight")
    ax[1].set_yscale('log')
    ax[1].set_xlabel(xaxis_label)
    ax[1].legend()
    ax[1].set_xlim(xlim)

    plt.tight_layout()
    plt.savefig(filename)
    plt.show()

    # Combined Data/MC ratio plot
    ratio_medium = contents_data[0] / contents_mc[0]
    ratio_tight = contents_data[1] / contents_mc[1]

    fig_ratio, ax_ratio = plt.subplots(figsize=(8, 6), dpi=150)
    ax_ratio.plot(x_centers, ratio_medium, 'o-', color='blue', markersize=4, label="Data/MC Medium")
    ax_ratio.plot(x_centers, ratio_tight, 'o-', color='green', markersize=4, label="Data/MC Tight")
    ax_ratio.axhline(1, color='black', linestyle='--')
    ax_ratio.set_ylabel("Data/MC Ratio")
    ax_ratio.set_xlabel(xaxis_label)
    ax_ratio.legend()
    ax_ratio.set_ylim((0, 2))
    ax_ratio.set_xlim(xlim)

    plt.tight_layout()
    plt.savefig(f"ratio_{filename}")
    plt.show()



# control_plots("d_tp_mc_minv", "d_tp_data_minv", "control_minv.png", "$M_{inv} [GeV]$", (0, 60))
# control_plots("d_tp_mc_aco", "d_tp_data_aco", "control_aco.png", "$Acoplanarity$", (0, 0.02))

def control_kin(histname, filename, xaxis_label, lab):
    plt.rcParams.update({'font.size': 15})
    # Open the ROOT file and extract the histogram
    file_uproot = uproot.open("output-basev2.root")
    hist = file_uproot[histname]
    content, edge = hist.to_numpy()
    y_err = hist.variances()
    y_err = np.sqrt(y_err)
    
    # Calculate bin centers and widths
    bin_centers = (edge[:-1] + edge[1:]) / 2
    bin_widths = (edge[1:] - edge[:-1]) / 2

    # Plot data as dots with x-errors
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    if edge[0] == 0:
        bin_centers=bin_centers[2:]
        content = content[2:]
        bin_widths = bin_widths[2:]
        y_err = y_err[2:]
    ax.errorbar(bin_centers, content, xerr=bin_widths, yerr=y_err, fmt='o', capsize=3, color='blue', label=lab)

    # Add labels and save the plot
    ax.set_xlabel(xaxis_label)

    ax.set_ylabel("Eff")
    ax.set_ylim(np.min(content)-0.1, 1.0)
    ax.legend(loc='lower right')
    plt.savefig(filename)
    plt.show()


# control_kin("d_eff_egamma_pt", "egamma_eff_pt.png", "$p_T^{truth}$ [GeV]", "rec Egam clu")
# control_kin("d_eff_pixtrk_pt", "pixtrk_eff_pt.png", "$p_T^{truth}$ [GeV]", "rec pixtrk")
# control_kin("d_eff_e_pt", "e_eff_pt.png", "$p_T^{truth}$ [GeV]", "rec electron")

# control_kin("d_eff_egamma_etaa", "egamma_eff_eta.png", "$\eta ^{truth}$", "rec Egam clu")
# control_kin("d_eff_pixtrk_eta", "pixtrk_eff_eta.png", "$\eta ^{truth}$", "rec pixtrk")
# control_kin("d_eff_e_eta", "e_eff_eta.png", "$\eta ^{truth}$", "rec electron")


def tagpro(hist_mc, hist_data, xlabel, save_file):
    file_uproot = uproot.open("output-basev2.root")
    mc = file_uproot[hist_mc]
    data = file_uproot[hist_data]
    mc_content, edges = mc.to_numpy()
    data_content, _ = data.to_numpy()
    data_content = data_content / np.sum(data_content)
    mc_content = mc_content / np.sum(mc_content)
    x = (edges[:-1] + edges[1:]) / 2
    x_err = (edges[1:] - edges[:-1]) / 2
    fig, ax = plt.subplots(figsize=(10,10), dpi=200)
    ax.errorbar(x, data_content, xerr=x_err, color="black", label="Data")
    ax.hist(edges[:-1], bins=edges, weights=mc_content, color='red', histtype='step', label="MC")
    ax.set_xlabel(xlabel)
    ax.legend(loc="best")
    ax.axvline(0.015, color='blue', linestyle='--')
    ax.set_ylabel("Events")
    ax.set_yscale("log")
    plt.savefig(save_file)
    plt.show()

#tagpro("d_tp_mc_aco", "d_tp_data_aco", "A(tag, probe)", "Aco.png")
# tagpro("d_tp_mc_minv", "d_tp_data_minv", "$M_{inv}^{tag,probe}$", "Minv.png")
# tagpro("d_dr_mc_probe", "d_dr_data_probe", "$\Delta R (rec,probe)$", "DR.png")


def comparison():
    file_uproot = uproot.open("output-basev2.root")
    hists_data_pt = ["d_rec_eff_data_pt_vl", "d_rec_eff_data_pt_lo", "d_rec_eff_data_pt_me", "d_rec_eff_data_pt_ti"]
    hists_mc_pt = ["d_rec_eff_mc_pt_vl", "d_rec_eff_mc_pt_lo", "d_rec_eff_mc_pt_me", "d_rec_eff_mc_pt_ti"]
    hists_data_eta = ["d_rec_eff_data_eta_vl", "d_rec_eff_data_eta_lo", "d_rec_eff_data_eta_me", "d_rec_eff_data_eta_ti"]
    hists_mc_eta = ["d_rec_eff_mc_eta_vl", "d_rec_eff_mc_eta_lo", "d_rec_eff_mc_eta_me", "d_rec_eff_mc_eta_ti"]
    bins_eta = [None]
    bins_pt = [None]
    
    for i in range(len(hists_data_eta)):
        hist1 = file_uproot[hists_data_pt[i]]
        hist2 = file_uproot[hists_mc_pt[i]]
        hist3 = file_uproot[hists_data_eta[i]]
        hist4 = file_uproot[hists_mc_eta[i]]
        var1, bins_pt = hist1.to_numpy()
        var2, _ = hist2.to_numpy()
        var3, bins_eta = hist3.to_numpy()
        var4, _ = hist4.to_numpy()
        hists_data_pt[i] = numpy.array(var1)
        hists_mc_pt[i] = numpy.array(var2)
        hists_data_eta[i] = numpy.array(var3)
        hists_mc_eta[i] = numpy.array(var4)
    
    bins_eta = numpy.array(bins_eta)
    bins_pt = numpy.array(bins_pt)
    colors = ["green", "red", "blue", "black"]

    # Separate legend elements
    op_legend_elements = [
        Line2D([0], [0], marker='o', color='green', label='VeryLoose', linestyle='None'),
        Line2D([0], [0], marker='o', color='red', label='Loose', linestyle='None'),
        Line2D([0], [0], marker='o', color='blue', label='Medium', linestyle='None'),
        Line2D([0], [0], marker='o', color='black', label='Tight', linestyle='None')
    ]
    data_mc_legend_elements = [
        Line2D([0], [0], marker='o', color='black', label='Data', markerfacecolor='black', linestyle='None'),
        Line2D([0], [0], marker='o', color='black', label='MC', markerfacecolor='none', markeredgecolor='black', linestyle='None')
    ]
    
    # Plot for pT
    fig, ax = plt.subplots(figsize=(8, 10), dpi=200)
    x = (bins_pt[:-1] + bins_pt[1:]) / 2
    x_err = (bins_pt[1:] - bins_pt[:-1]) / 2
    for i in range(4):
        ax.errorbar(x=x, y=hists_data_pt[i], xerr=x_err, fmt='o', ecolor=colors[i], markerfacecolor=colors[i], markeredgecolor=colors[i])
        ax.errorbar(x=x, y=hists_mc_pt[i], xerr=x_err, fmt='o', ecolor=colors[i], markerfacecolor='none', markeredgecolor=colors[i])
    
    ax.set_xlabel("$p_T^{probe} [GeV]$")
    ax.set_ylabel("Identification efficiency")
    ax.set_ylim((0, 1))
    
    # Adding legends side by side
    op_legend = ax.legend(handles=op_legend_elements, loc='lower right', bbox_to_anchor=(1, 0), title="Operating Points")
    ax.add_artist(op_legend)  # Add the first legend
    ax.legend(handles=data_mc_legend_elements, loc='lower right', bbox_to_anchor=(0.75, 0), title="Data/MC")
    
    plt.savefig("OP_comp_pt.png")
    plt.show()

    # Plot for eta
    fig, ax = plt.subplots(figsize=(8, 10), dpi=200)
    x = (bins_eta[:-1] + bins_eta[1:]) / 2
    x_err = (bins_eta[1:] - bins_eta[:-1]) / 2
    for i in range(4):
        ax.errorbar(x=x, y=hists_data_eta[i], xerr=x_err, fmt='o', ecolor=colors[i], markerfacecolor=colors[i], markeredgecolor=colors[i])
        ax.errorbar(x=x, y=hists_mc_eta[i], xerr=x_err, fmt='o', ecolor=colors[i], markerfacecolor='none', markeredgecolor=colors[i])
    
    ax.set_xlabel("$\eta^{probe}$")
    ax.set_ylabel("Identification efficiency")
    ax.set_ylim((0, 1))
    
    # Adding legends side by side
    op_legend = ax.legend(handles=op_legend_elements, loc='lower right', bbox_to_anchor=(1, 0), title="Operating Points")
    ax.add_artist(op_legend)  # Add the first legend
    ax.legend(handles=data_mc_legend_elements, loc='lower right', bbox_to_anchor=(0.75, 0), title="Data/MC")
    
    plt.savefig("OP_comp_eta.png")
    plt.show()


#comparison()

# file = uproot.open("output-basev2.root")
# hist_data = file["d_tp_data_pt"]
# hist_mc = file["d_tp_mc_pt"]
# content_mc, edges = hist_mc.to_numpy()
# content_data, _ = hist_data.to_numpy()
# plt.figure(dpi=200)
# x_centers = (edges[:-1] + edges[1:]) / 2
# content_mc = numpy.array(content_mc)
# content_data = numpy.array(content_data)
# content_mc = content_mc / numpy.sum(content_mc)
# content_data = content_data / numpy.sum(content_data)
# edges = numpy.array(edges)
# plt.hist(edges[:-1], bins=edges, weights=content_mc, color='red', histtype='step', label="MC")
# plt.scatter(x_centers, content_data, color="black", label="Data", s=10.)
# plt.legend(loc="best")
# plt.xlabel("$p_T [GeV]$")
# plt.yscale("log")
# plt.savefig("d0.png")
# plt.show()


# file_tig = uproot.open("output-tigv2.root")
# file_base = uproot.open("output-basev2.root")
# hist_tig = file_tig["d_ratio_eta"]
# hist_base = file_base["d_ratio_eta"]
# content_tig, edges = hist_tig.to_numpy()
# content_base, _ = hist_base.to_numpy()
# content_tig = numpy.array(content_tig)
# content_base = numpy.array(content_base)
# edges = numpy.array(edges)
# plt.figure(dpi=200)
# plt.hist(edges[:-1], bins=edges, weights=content_base, color='red', histtype='step', label="Medium")
# plt.hist(edges[:-1], bins=edges, weights=content_tig, color='blue', histtype='step', label="Tight")
# plt.legend(loc="best")
# plt.xlabel("$\eta$")
# plt.ylabel("SF")
# plt.ylim(0.8, 2.0)
# plt.savefig("SF_comp.png")

