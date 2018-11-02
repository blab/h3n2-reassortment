import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from parse_epi_data import parse_ili, parse_subtype_distribution
from collections import defaultdict
import sys

from matplotlib.colors import ListedColormap
from matplotlib.contour import ContourSet
import matplotlib.cm as cm

figure_output_dir = "../figures/flu_epi/"

plt.ion()
sns.palplot(sns.color_palette("Blues"))

sns.set_style('ticks')
h1 = 'A (H1)'
h3 = 'A (H3)'
h1pdm = 'A (2009 H1N1)'
A_not = 'A (Subtyping not Performed)'
B = 'B'
subtypes = [h1, h1pdm, h3, B]
cols = {k:c for c, k in zip(sns.color_palette(n_colors=6), subtypes)}
fs=16

ili, ili_seasons = parse_ili()
st, dominant_strain = parse_subtype_distribution()

season_colors = sns.color_palette("Reds", n_colors=len(ili_seasons))
ordered_seasons = sorted([ k for k in ili_seasons.keys()])
###
# PLOT ILI FOR ALL YEARS COLORED BY YEAR
###

Z = [[0,0],[0,0]]
levels = range(len(ordered_seasons)+1)
CS3 = plt.contourf(Z, levels, cmap=ListedColormap(season_colors.as_hex()))
plt.clf()

fig = plt.figure(figsize=(11,8))
lines = [ plt.plot(ili_seasons[ordered_seasons[i]][:,0],
            ili_seasons[ordered_seasons[i]][:,1],
            c=season_colors[i], ls='-',
            lw = 4 if ordered_seasons[i]=="2017-2018" else 2,
            alpha= 1.0 if ordered_seasons[i]=="2017-2018" else 0.5) for i in range(len(ili_seasons)) ]

cbar = plt.colorbar(CS3)
cbar.ax.get_yaxis().set_ticks([])
for j, lab in enumerate(ordered_seasons):
    cbar.ax.text(2.9, (j / 21.)+.02, lab, ha='center', va='center', fontsize=fs, fontweight='bold')

ax = plt.gca()
ax.text(-8.3,7.2, "2009-2010\n(H1N1pdm)", fontsize=fs)
ax.text(7,7.2, "2017-2018\n(H3N2)", fontsize=fs)
# plt.legend(fontsize=fs)
# plt.legend.remove()
for tick in ax.get_xticklabels():
    tick.set_fontsize(1.75*fs)
    # tick.set_fontweight('bold')
for tick in ax.get_yticklabels():
    tick.set_fontsize(1.75*fs)
    # tick.set_fontweight('bold')
plt.tick_params(labelsize=0.8*fs)
plt.ylabel('ILI [%]', fontsize=1.25*fs, fontweight='bold')
# plt.xlim(-10, 40) for full year view
plt.xlim(-10, 20)
# plt.xticks([-10, 0, 10, 20, 30, 40], map(str,[42, 0, 10, 20, 30, 40] ))
plt.xticks([-10, 0, 10, 20], map(str,[42, 0, 10, 20] ))
plt.xlabel("week of year", fontsize=1.25*fs, fontweight='bold')
# plt.savefig('{}ili_1997_2018.png'.format(figure_output_dir))
plt.savefig('{}ili_1997_2018.svg'.format(figure_output_dir), format="svg")
sys.exit()

##############################################
###                                        ###
###          ALL FOLLOWING UNUSED          ###
###                                        ###
##############################################

###
# subtypes through time
###
plt.figure(figsize=(12,6))
for ii,tri in enumerate(subtypes):
    plt.plot_date(st.date, st[tri], ls='-', label = tri, c=cols[tri])

plt.legend(loc=2, fontsize=fs, ncol=2)
plt.ylabel('# of positive tests per subtype', fontsize=fs)
plt.tick_params(labelsize=0.8*fs)
plt.savefig('flu_subtypes_positive_tests_linear.png'.format(figure_output_dir))
plt.ylim(5,20000)
plt.yscale('log')
plt.savefig('flu_subtypes_positive_tests.png'.format(figure_output_dir))


####
# fit a Gaussian to the peak and extract peak width, position
# and the preseason growth rate
# ##
from scipy.optimize import minimize

def peak_fit_func(params, x, y):
    res = (y - params[0]*np.exp(-0.5*(x-params[1])**2/params[2]**2))**2
    return np.sum(res)

def slope_fit_func(params, x, y):
    #res = (np.log(y+1) - np.log(max(1,params[0]))-params[1]*x)**2*np.maximum(0,np.log(y+1))
    res = (y - params[0]*np.exp(params[1]*x))**2/np.maximum(1,y)
    return np.sum(res)

slopes_weighted = defaultdict(dict)
rev_slopes = defaultdict(dict)
for ii,tri in enumerate(subtypes):
    for season in np.unique(st["season"]):
        if season=="2009-2010" or (tri=='A (2009 H1N1)' and season=='2008-2009'):
            continue # ignore pandemic season
        data = st.loc[st["season"]==season,(tri,"season_week")]
        if data[tri].max()>400:
            shift = data["season_week"][data[tri].argmax()]
            val = data[tri].max()
            sol = minimize(peak_fit_func, [val, shift, 4],
                           args=(data["season_week"], data[tri]))
            val, shift, width = sol['x']
            if width<0:
                print(tri, season, "has negative width")
                continue

            upper = (data[tri]>val*0.5).argmax()
            lower = (data[tri]>70).argmax()
            dslice = data.loc[lower:upper]
            # generate starting values by linear regresson
            res = linregress(dslice['season_week']-shift, np.log(dslice[tri]))
            # fit non-linear model that weights points by number of observations
            sol = minimize(slope_fit_func, [np.exp(res.intercept), res.slope],
                           args=(dslice["season_week"]-shift, dslice[tri]))
            slopes_weighted[tri][season] = [sol["x"][1], val, shift, width]

            past_peak = data.loc[data[tri].argmax():]
            lower = (past_peak[tri]<val*0.5).argmax()
            upper = (past_peak[tri]<100).argmax()
            dslice = past_peak.loc[lower:upper]
            sol = minimize(slope_fit_func, [val, -0.2],
                           args=(dslice["season_week"]-shift, dslice[tri]))
            rev_slopes[tri][season] = [sol["x"][1], val, shift, width]


slopes_weighted_array = {}
for tri in slopes_weighted:
    slopes_weighted_array[tri]=np.array(list(slopes_weighted[tri].values()))
rev_slopes_array = {}
for tri in rev_slopes:
    rev_slopes_array[tri]=np.array(list(rev_slopes[tri].values()))

####
# plot normalized number of positive tests and preseason growth rates
###
#fig, axs = plt.subplots(1,3, figsize=(15,6))
fig, axs = plt.subplots(1,2, figsize=(12,6))
already_labeled = {k:False for k in subtypes}
for ii,tri in enumerate(subtypes):
    for season in np.unique(st["season"]):
        data = st.loc[st["season"]==season,(tri,"season_week")]
        if season in slopes_weighted[tri]:
            slope, val, shift, width = slopes_weighted[tri][season]
            axs[0].plot(data["season_week"]-shift, data[tri]/val,
                        alpha=0.5, c=cols[tri], ls='-', lw=3,
                        label = "" if already_labeled[tri] else tri)
            already_labeled[tri]=True

axs[0].set_ylabel('normalized positive tests',size=fs)
axs[0].set_xlabel('week relative to peak',size=fs)
axs[0].set_ylim(0.001,2)
axs[0].set_xlim(-10,20)
axs[0].set_yscale('log')
axs[0].legend()
axs[0].tick_params(labelsize=0.8*fs)

sns.boxplot(data=[slopes_weighted_array[k][:,0] for k in subtypes], ax=axs[1], color="0.9")
#            color = [cols[k] for k in subtypes])
sns.stripplot(data=[slopes_weighted_array[k][:,0] for k in subtypes],
                    color=".4", size=8, ax=axs[1])

for tri in subtypes:
    if "2017-2018" in slopes_weighted[tri]:
        axs[1].scatter([subtypes.index(tri)], [slopes_weighted[tri]["2017-2018"][0]], c='r', marker='+', s=200)
        print([subtypes.index(tri)], [slopes_weighted[tri]["2017-2018"][0]])

axs[1].set_ylim([0,0.8])
axs[1].set_ylabel('growth rate [1/week]', size=fs)
axs[1].set_xticklabels(subtypes, rotation=30,size=fs*0.8)
axs[1].tick_params(labelsize=0.8*fs)


# sns.boxplot(data=[rev_slopes_array[k][:,0] for k in subtypes], ax=axs[2], color="0.9")
# #            color = [cols[k] for k in subtypes])
# sns.stripplot(data=[rev_slopes_array[k][:,0] for k in subtypes],
#                     color=".4", size=8, ax=axs[2])

# axs[2].set_ylim([-0.6,0])
# axs[2].set_ylabel('shrink rate [1/week]', size=fs)
# axs[2].set_xticklabels(subtypes, rotation=30,size=fs*0.8)
# axs[2].tick_params(labelsize=0.8*fs)

plt.tight_layout()
plt.savefig('{}flu_growth_rate.png'.format(figure_output_dir))

###
# plot peak week and width
###
fig, axs = plt.subplots(1,2, figsize=(12,6))

sns.boxplot(data=[slopes_weighted_array[k][:,3] for k in subtypes], ax=axs[0], color="0.9")
#            color = [cols[k] for k in subtypes])
sns.stripplot(data=[slopes_weighted_array[k][:,3] for k in subtypes],
                    color=".4", size=8, ax=axs[0])

sns.boxplot(data=[slopes_weighted_array[k][:,2] for k in subtypes], ax=axs[1], color="0.9")
#            color = [cols[k] for k in subtypes])
sns.stripplot(data=[slopes_weighted_array[k][:,2] for k in subtypes],
                    color=".4", size=8, ax=axs[1])

axs[1].set_ylim((-5,20))
axs[1].set_ylabel('peak week', size=fs)
axs[1].set_xticklabels(subtypes, rotation=30,size=fs*0.8)
axs[1].tick_params(labelsize=0.8*fs)

axs[0].set_ylim((0,7))
axs[0].set_ylabel('peak standard dev [weeks]', size=fs)
axs[0].set_xticklabels(subtypes, rotation=30,size=fs*0.8)
axs[0].tick_params(labelsize=0.8*fs)

plt.tight_layout()
plt.savefig('{}flu_peak_properties.png'.format(figure_output_dir))
