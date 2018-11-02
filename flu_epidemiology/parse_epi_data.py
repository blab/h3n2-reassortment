import pandas as pd
import numpy as np
from collections import defaultdict

ili_fname = '../data/flu_epi_data/ILINet_new.csv'
subtype_fname_prior_2016 = '../data/flu_epi_data/WHO_NREVSS_Combined_prior_to_2015_16.csv'
subtype_fname_since_2016_PH = '../data/flu_epi_data/WHO_NREVSS_Public_Health_Labs.csv'

def year_week_to_date(year, week):
    import datetime
    return datetime.date.fromordinal(datetime.date(year, 1,1).toordinal()+(week-1)*7)


def year_week_to_season(year, week):
	if week<40:
		return "%d-%d"%(year-1, year), week
	else:
		return "%d-%d"%(year, year+1), week-52


def to_float(x):
    if x=='X':
        return np.nan
    else:
        try:
            return float(x)
        except:
            return np.nan


def parse_ili():
	ili = pd.read_csv(ili_fname, skiprows=1)
	ili_seasons=defaultdict(list)
	for ri,row in ili.iterrows():
		season, sweek = year_week_to_season(row.YEAR, row.WEEK)
		ili_seasons[season].append((sweek, row["%UNWEIGHTED ILI"]))
	for season in ili_seasons:
		ili_seasons[season] = np.array(ili_seasons[season])
	return ili, ili_seasons


def parse_subtype_distribution():
	h1 = 'A (H1)'
	h3 = 'A (H3)'
	h1pdm = 'A (2009 H1N1)'
	A_not = 'A (Subtyping not Performed)'
	B = 'B'

	st_prior_2016 = pd.read_csv(subtype_fname_prior_2016, skiprows=1)
	st_past_2016 = pd.read_csv(subtype_fname_since_2016_PH, skiprows=1)

	st_combined = pd.concat((st_prior_2016, st_past_2016))
	st_combined = st_combined.replace(np.nan, 0.0)
	st_combined["season"] = pd.Series([year_week_to_season(row.YEAR, row.WEEK)[0]
								for ri,row in st_combined.iterrows()], index=st_combined.index)
	st_combined["season_week"] = pd.Series([year_week_to_season(row.YEAR, row.WEEK)[1]
								for ri,row in st_combined.iterrows()], index=st_combined.index)
	st_combined["date"] = pd.Series([year_week_to_date(row.YEAR, row.WEEK)
								for ri,row in st_combined.iterrows()], index=st_combined.index)

	missing = {}
	for strain in [h3, h1, h1pdm]:
		missing[strain] = st_combined[A_not]*st_combined[strain]/(st_combined[h3]+st_combined[h1]+st_combined[h1pdm]+0.01)
	for strain in missing:
		st_combined[strain]+=missing[strain]

	dominant_strain = {}
	for season in np.unique(st_combined["season"]):
		dominant_strain[season] = st_combined.loc[st_combined["season"]==season,(h1,h1pdm,h3,B)].sum(axis=0).argmax()

	return st_combined, dominant_strain
