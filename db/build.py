#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
This script builds the AstroRV catalogue.

This script is designed to be ran on ipython
"""
print("Building AstroRV catalogue")

# =============================================================================
# EXTERNAL PACKAGES
# =============================================================================
import sys,os,urllib,collections
sys.path.append("/home/local-python/lib/python3.5/site-packages")
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from sys import argv
exit=sys.exit

if len(argv)==1:qread=1
else:qread=int(argv[1])

# =============================================================================
# DIRECTORIES
# =============================================================================
#COMMON DIR
SRC_DIR="src/"

#GAIA SOURCES
TGAS_DIR="src/Astro/"

#SIMBAD SOURCES
SIMBAD_DIR="src/Astro/"

#TYCHO2/HIPPARCOS SOURCES
HIPTYC_DIR="src/Astro/"

#RADIAL VELOCITY SOURCES
RV_DIR="src/RV/"

if qread:
    # =============================================================================
    # DOWNLOAD GAIA
    # =============================================================================
    print("Downloading data from Gaia..")
    url_dir = "http://cdn.gea.esac.esa.int/Gaia/tgas_source/csv/"
    for i in range(16):
        file = "TgasSource_000-000-0" + str(i).zfill(2) + ".csv.gz"
        url = url_dir + file
        filename = TGAS_DIR + file
        if not os.path.isfile(filename):
            print("\tDownloading", file)
            urllib.request.urlretrieve(url, filename)
        else:
            print("\tFile ",file," already downloaded")

    # =============================================================================
    # READ GAIA DATA
    # =============================================================================
    print("Reading and building Gaia database...")
    cols_gaia = ["hip", "tycho2_id", "ra", "ra_error", "dec", "dec_error", "parallax", "parallax_error", "pmra", "pmra_error",         "pmdec", 
                 "pmdec_error", "ra_dec_corr", "ra_parallax_corr", "ra_pmra_corr", "ra_pmdec_corr", "dec_parallax_corr",         "dec_pmra_corr", 
                 "dec_pmdec_corr", "parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr",         "phot_g_mean_flux", "phot_g_mean_flux_error", 
                 "phot_g_mean_mag", "l", "b", "ecl_lon", "ecl_lat"]

    for i in range(16):
        filename = TGAS_DIR + "TgasSource_000-000-0" + str(i).zfill(2) + ".csv.gz"
        if i == 0:
            print("\tReading", filename)
            gaia = pd.read_csv(filename, usecols=cols_gaia)
        else:
            print("\tReading", filename)
            DRx = pd.read_csv(filename, usecols=cols_gaia)
            gaia = gaia.append(DRx)

    gaia = pd.DataFrame(gaia)
    gaia_hip = gaia[gaia.hip.notnull()]
    gaia_tyc = gaia[gaia.tycho2_id.notnull()]

    print("\tGaia: Subset Hipparcos:", len(gaia_hip))
    print("\tGaia: Subset Tycho-2:", len(gaia_tyc))
    print("\tTotal Gaia:", len(gaia))

    # =============================================================================
    # READ HIPPARCOS DATA
    # =============================================================================
    print("Reading and building Hipparcos database...")
    # Información disponible en: http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=I%2F239&target=readme&#sRM2.1
    names_hip={}
    names_hip[1] = "hip"
    names_hip[8] = "ra_Hip"
    names_hip[9] = "dec_Hip"
    names_hip[11] = "parallax_Hip"
    names_hip[12] = "pmra_Hip"
    names_hip[13] = "pmdec_Hip"
    names_hip[14] = "ra_error_Hip"
    names_hip[15] = "dec_error_Hip"
    names_hip[16] = "parallax_error_Hip"
    names_hip[17] = "pmra_error_Hip"
    names_hip[18] = "pmdec_error_Hip"
    names_hip[19] = "ra_dec_corr_Hip"
    names_hip[20] = "ra_parallax_corr_Hip"
    names_hip[21] = "dec_parallax_corr_Hip"
    names_hip[22] = "ra_pmra_corr_Hip"
    names_hip[23] = "dec_pmra_corr_Hip"
    names_hip[24] = "parallax_pmra_corr_Hip"
    names_hip[25] = "ra_pmdec_corr_Hip"
    names_hip[26] = "dec_pmdec_corr_Hip"
    names_hip[27] = "parallax_pmdec_corr_Hip"
    names_hip[28] = "pmra_pmdec_corr_Hip"
    names_hip[71] = "HenryDraperId_Hip"
    names_hip[5] = "Vmag_Hip"
    #names_hip[76] = "sptype_Hip"

    names_hip = collections.OrderedDict(sorted(names_hip.items()))

    hipparcos = pd.read_csv(HIPTYC_DIR+"hip_main.dat.gz", delimiter="|", usecols = names_hip.keys(),  names = names_hip.values())
    hipparcos = pd.DataFrame(hipparcos)

    n1 = len(hipparcos)
    hipparcos["ra_Hip"] = pd.to_numeric(hipparcos["ra_Hip"], errors="coerce")
    hipparcos["dec_Hip"] = pd.to_numeric(hipparcos["dec_Hip"], errors="coerce")
    hipparcos["parallax_Hip"] = pd.to_numeric(hipparcos["parallax_Hip"], errors="coerce")
    hipparcos["Vmag_Hip"] = pd.to_numeric(hipparcos["Vmag_Hip"], errors="coerce")
    hipparcos.dropna(subset=["ra_Hip", "dec_Hip", "parallax_Hip"], inplace=True)
    hipparcos["source"]="hipparcos"

    n2 = len(hipparcos)
    print("\tObjects discarded:", n1-n2)
    print("\tFinal size of database:",n2)

    # =============================================================================
    # READ TYCHO DATA
    # =============================================================================
    print("Reading and building Tycho database...")
    # Información disponible en: http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=I%2F239&target=readme&#sRM2.13
    names_tyc={}
    names_tyc[1] = "tycho2_id"
    names_tyc[8] = "ra_Tyc"
    names_tyc[9] = "dec_Tyc"
    names_tyc[11] = "parallax_Tyc"
    names_tyc[12] = "pmra_Tyc"
    names_tyc[13] = "pmdec_Tyc"
    names_tyc[14] = "ra_error_Tyc"
    names_tyc[15] = "dec_error_Tyc"
    names_tyc[16] = "parallax_error_Tyc"
    names_tyc[17] = "pmra_error_Tyc"
    names_tyc[18] = "pmdec_error_Tyc"
    names_tyc[19] = "ra_dec_corr_Tyc"
    names_tyc[20] = "ra_parallax_corr_Tyc"
    names_tyc[21] = "dec_parallax_corr_Tyc"
    names_tyc[22] = "ra_pmra_corr_Tyc"
    names_tyc[23] = "dec_pmra_corr_Tyc"
    names_tyc[24] = "parallax_pmra_corr_Tyc"
    names_tyc[25] = "ra_pmdec_corr_Tyc"
    names_tyc[26] = "dec_pmdec_corr_Tyc"
    names_tyc[27] = "parallax_pmdec_corr_Tyc"
    names_tyc[28] = "pmra_pmdec_corr_Tyc"
    names_tyc[53] = "HenryDraperId_Tyc"
    names_tyc[5] = "Vmag_Tyc"

    names_tyc = collections.OrderedDict(sorted(names_tyc.items()))

    tycho = pd.read_csv(HIPTYC_DIR+"tyc_main.dat", delimiter="|", usecols = names_tyc.keys(), names = names_tyc.values())
    tycho = pd.DataFrame(tycho)

    tycho["a"], tycho["b"], tycho["c"] = tycho["tycho2_id"].str.split().str
    tycho["tycho2_id"] = tycho["a"] + "-" + tycho["b"] + "-" + tycho["c"]
    del tycho["a"], tycho["b"], tycho["c"]

    n1 = len(tycho)
    tycho["ra_Tyc"] = pd.to_numeric(tycho["ra_Tyc"], errors="coerce")
    tycho["dec_Tyc"] = pd.to_numeric(tycho["dec_Tyc"], errors="coerce")
    tycho["parallax_Tyc"] = pd.to_numeric(tycho["parallax_Tyc"], errors="coerce")
    tycho["Vmag_Tyc"] = pd.to_numeric(tycho["Vmag_Tyc"], errors="coerce")
    tycho.dropna(subset=["ra_Tyc", "dec_Tyc", "parallax_Tyc"], inplace=True)
    tycho["source"]="tycho"

    n2 = len(tycho)
    print("\tObjects discarded:", n1-n2)
    print("\tFinal size of database:",n2)

    # =============================================================================
    # READ SIMBAD DATA
    # =============================================================================
    print("Reading and building Simbad database...")

    cols = ["typedident","identifier", "radvel", "coord1(ICRS,J2000/2000)", "plx", "pm", "MagV", "spec.type"]
    simbad = pd.read_csv(SIMBAD_DIR+"simbad.csv", usecols=cols, delimiter="|")
    simbad = pd.DataFrame(simbad)

    simbad["hip"] = simbad["typedident"].map(lambda x: str(x)[4:]).astype(float)
    del simbad["typedident"]

    simbad["plx"] = pd.to_numeric(simbad["plx"], errors="coerce")   

    simbad["coord1(ICRS,J2000/2000)"] = simbad["coord1(ICRS,J2000/2000)"].str.strip()
    simbad["ra_h"], simbad["ra_m"], simbad["ra_s"], simbad["dec_h"], simbad["dec_m"], simbad["dec_s"] =     simbad["coord1(ICRS,J2000/2000)"].str.split(" ").str
    simbad["ra_simbad"] = simbad["ra_h"].astype(float)*15 + simbad["ra_m"].astype(float)/60 + simbad["ra_s"].astype(float)/3600
    simbad["dec_simbad"] = np.sign(simbad["dec_h"].astype(float)) * (     np.abs(simbad["dec_h"].astype(float)) + simbad["dec_m"].astype(float)/60 + simbad["dec_s"].astype(float)/3600 )

    del simbad["coord1(ICRS,J2000/2000)"]
    del simbad["ra_h"], simbad["ra_m"], simbad["ra_s"], simbad["dec_h"], simbad["dec_m"], simbad["dec_s"]

    simbad["pm"] = simbad["pm"].str.strip()
    simbad["pmra_simbad"], simbad["pmdec_simbad"] = simbad["pm"].str.split(" ").str
    del simbad["pm"]

    n1 = len(simbad)
    simbad.dropna(subset=["ra_simbad", "dec_simbad", "plx"], inplace=True)      
    n2 = len(simbad)
    print("\tObjects discarded:", n1-n2)

    simbad["radvel"] = simbad["radvel"].str.strip()
    simbad["radvel"] = simbad["radvel"].replace("~", np.nan)

    simbad = simbad.rename(columns={"identifier": "name_simbad"})
    simbad = simbad.rename(columns={"plx": "parallax_simbad"})
    simbad = simbad.rename(columns={"spec.type": "sptype_simbad"})
    simbad = simbad.rename(columns={"radvel": "radial_vel_simbad"})
    simbad = simbad.rename(columns={"MagV": "Vmag_simbad"})
    simbad["source"]="simbad"

    print("\tFinal size of database:",len(simbad))

def buildCat(gaia,hipparcos,tycho,simbad):
    # =============================================================================
    # MERGING ASTROMETRIC DATABASES
    # =============================================================================
    print("Merging astrometric databases...")
    gaia_hip = gaia[gaia.hip.notnull()]
    gaia_tyc = gaia[gaia.tycho2_id.notnull()]

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #COMPACT TABLE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #"""
    print("\tCompact table...")
    gaia["source"] = "gaia"
    exclusive_hip = hipparcos[~hipparcos.hip.isin(gaia_hip.hip)]
    exclusive_hip.rename(columns=lambda x: x.replace("_Hip", ""), inplace=True)
    exclusive_hip["source"] = "hipparcos"
    db = pd.concat([gaia, exclusive_hip], axis=0)
    print("\tHipparcos contribution:",len(exclusive_hip))

    exclusive_tyc = tycho[~tycho.tycho2_id.isin(gaia_tyc.tycho2_id)]
    exclusive_tyc.rename(columns=lambda x: x.replace("_Tyc", ""), inplace=True)
    exclusive_tyc["source"] = "tycho"
    db = pd.concat([db, exclusive_tyc], axis=0)
    print("\tTycho contribution:", len(exclusive_tyc))

    exclusive_simbad = simbad[~simbad.hip.isin(db.hip)]
    exclusive_simbad.rename(columns=lambda x: x.replace("_simbad", ""), inplace=True)
    exclusive_simbad = exclusive_simbad.loc[:,["hip", "ra", "dec", "parallax", "pmra", "pmdec"]]
    exclusive_simbad["source"] = "simbad"
    db = pd.concat([db, exclusive_simbad], axis=0)
    print("\tSimbad contribution:", len(exclusive_simbad))

    db = db.merge(simbad.loc[:,["hip", "name_simbad","sptype_simbad","radial_vel_simbad","Vmag_simbad"]],               left_on="hip", right_on="hip", how="left")

    cols_gaia=list(gaia.columns)
    cols = ["source"] + cols_gaia + ["Vmag","HenryDraperId"] + ["name_simbad","sptype_simbad","radial_vel_simbad","Vmag_simbad"]

    db = db[cols]
    db = db.reset_index(drop=True)

    print("\tFinal size:", len(db))
    
    dbfile=TGAS_DIR+"AstroComp.csv"
    print("\tSaving ",dbfile)
    db.to_csv(dbfile,index=False)
    print("\tDone.")

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #FULL TABLE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print("\tFull table...")

    n1 = len(gaia)
    print("\tInitial size:", n1)

    database = gaia.merge(hipparcos, left_on="hip", right_on="hip", how="outer")
    n2 = len(database)
    print("\tHipparcos size:", len(hipparcos), ". Final size of database:", n2, ". Contribution:", n2-n1)

    database = database.merge(tycho, left_on="tycho2_id", right_on="tycho2_id", how="outer")
    n3 = len(database)
    print("\tTycho size:", len(tycho), ". Final size of database:", n3, ". Contribution:", n3-n2)

    database = database.merge(simbad, left_on="hip", right_on="hip", how="outer")
    n4 = len(database)
    print("\tSimbad size:", len(simbad), ". Final size of database:", n4, ". Contribution:", n4-n3)

    sources=database[['source_x','source_y']].fillna('')
    #sources.columns=['s1','s2','s3','s4']
    #database["source"]=sources["s1"].astype(str)+":"+sources["s2"]+":"+sources["s3"]+":"+sources["s4"]
    sources.columns=['s1','s2']
    database["source"]=sources["s1"].astype(str)+":"+sources["s2"]
    del(database["source_x"])
    del(database["source_y"])

    database = database.reset_index(drop=True)

    dbfile=TGAS_DIR+"Astro.csv"
    print("\tSaving ",dbfile,)
    database.to_csv(dbfile,index=False)
    print("\tDone.")

    # =============================================================================
    # RADIAL VELOCITY DATABASES
    # =============================================================================
    data=dict()
    match=dict()
    RVcat=pd.DataFrame()
    srcdir=RV_DIR

    ###################################################
    #READ CATALOGUES
    ###################################################
    #MALDONADO2010
    #J/A+A/521/A12/table1
    name="Maldonado2010.tsv"
    print("\tBuilding catalogue %s..."%name)
    comments=list(range(82))+[83,84]
    data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
    cats=['HIP']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="_RAJ2000",DEJ2000="_DEJ2000")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["RV"]
    df["eRV"]=data[name]["e_RV"]
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    #WEB1995
    #III/213/catalog
    name="Web1995-HIP.csv"
    print("\tBuilding catalogue %s..."%name)
    data[name]=pd.read_csv(srcdir+name)
    cats=['HIP']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="_RAJ2000",DEJ2000="_DEJ2000")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["RV"].map(str)
    df["eRV"]=data[name]["e_RV"].map(str)
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    #WEB1995-TYC2
    #III/213/catalog
    name="Web1995-TYC2.csv"
    print("\tBuilding catalogue %s..."%name)
    data[name]=pd.read_csv(srcdir+name)
    cats=['TYC2','HIP']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="_RAJ2000",DEJ2000="_DEJ2000")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["RV"].map(str)
    df["eRV"]=data[name]["e_RV"].map(str)
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    #GCS2011
    #J/A+A/530/A138/catalog
    name="GCS2011.tsv"
    print("\tBuilding catalogue %s..."%name)
    comments=list(range(174))+[175,176]
    data[name]=pd.read_csv(srcdir+name,sep="|",skiprows=comments)
    cats=['HIP']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="_RAJ2000",DEJ2000="_DEJ2000")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["RV"].map(str)
    df["eRV"]=data[name]["e_RV"].map(str)
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    #RAVE-DR5
    #III/279/rave_dr5
    name="RAVE-DR5.tsv"
    print("\tBuilding catalogue %s..."%name)
    comments=list(range(78))+[79,80]
    data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
    cats=['TYCHO2']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    data[name]["TYC2"]=data[name]["TYCHO2"]
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="RAJ2000",DEJ2000="DEJ2000")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["HRV"].map(str)
    df["eRV"]=data[name]["e_HRV"].map(str)
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    #PULKOVO
    #III/252/table8
    name="Pulkovo.tsv"
    print("\tBuilding catalogue %s..."%name)
    comments=list(range(61))+[62,63]
    data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
    cats=['HIP']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="_RA",DEJ2000="_DE")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["RV"].map(str)
    df["eRV"]=data[name]["eRV"].map(str)
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    #FAMAEY2005
    #J/A+A/430/165/tablea1
    name="Famaey2005.tsv"
    print("\tBuilding catalogue %s..."%name)
    comments=list(range(118))+[119,120]
    data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
    cats=['HIP']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="_RAJ2000",DEJ2000="_DEJ2000")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["RV"].map(str)
    df["eRV"]=data[name]["e_RV"].map(str)
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    #BB2000
    #III/213/catalogue
    name="BB2000.csv"
    print("\tBuilding catalogue %s..."%name)
    data[name]=pd.read_csv(srcdir+name)
    cats=['TYC2','HIP']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="_RAJ2000",DEJ2000="_DEJ2000")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["RV"].map(str)
    df["eRV"]=data[name]["e_RV"].map(str)
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    #MALARODA2012
    #III/249/catalog
    name="Malaroda2012.csv"
    print("\tBuilding catalogue %s..."%name)
    data[name]=pd.read_csv(srcdir+name)
    cats=['TYC2','HIP']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="_RAJ2000",DEJ2000="_DEJ2000")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["RV"].map(str)
    df["eRV"]=1.0;df["eRV"]=df["eRV"].map(str) #TYPICAL VALUE FOR OTHER CATALOGUES
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    #GALAH I
    #J/MNRAS/465/3203/catal
    name="Galah.tsv"
    print("\tBuilding catalogue %s..."%name)
    comments=list(range(54))+[55,56]
    data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
    cats=['TYC2']
    for cname in cats:
        data[name][cname]=data[name][cname].fillna('')
        data[name][cname]=data[name][cname].map(str)
    dfstr=data[name].select_dtypes(['object'])
    data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    print("\tNumber of objects in %s:"%name,len(data[name]))
    for cat in cats:
        cond=data[name][cat]!=''
        print("\tNumber of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
    #STORING RESULTS
    df=pd.DataFrame()
    if not 'TYC2' in data[name].columns:df["TYC2"]=''
    else:df["TYC2"]=data[name]["TYC2"]
    if not 'HIP' in data[name].columns:df["HIP"]=''
    else:df["HIP"]=data[name]["HIP"].apply(lambda x:x.replace('.0',''))
    COORDS=dict(RAJ2000="RAJ2000",DEJ2000="DEJ2000")
    for C in COORDS.keys():df[C]=data[name][COORDS[C]]
    df["RV"]=data[name]["RV"].map(str)
    df["eRV"]=0.6;df["eRV"]=df["eRV"].map(str) #Martell et al. (2017)
    df["CAT"]=name
    #ZERO ERRORS
    cond=df.RV==''
    df=df.drop(df.index[cond])
    df.RV=df.RV.map(float)
    df.eRV[df.eRV=='']=0.0
    df.eRV=df.eRV.map(float)
    med=df.eRV[df.eRV>0].median()
    print("\tMedian error:",med)
    print("\tNumber of entries with zero error:",len(df.eRV[df.eRV==0]))
    df.eRV[df.eRV==0]=med
    #RESULTING SIZE
    print("\tFiltered catalogue:",len(df))
    #FILLNA
    RVcat=RVcat.append(df.fillna(''))

    ###################################################
    #COMPILING FULL TABLE
    ###################################################
    print("\tCompiling final catalogue...")
    RVcat=RVcat.rename(columns={'TYC2':'tycho2_id','HIP':'hip'})
    print("\tCompiling final catalogue...")
    RVcatf=RVcat.reset_index(drop=True)
    print("\tNumber of unfiltered entries:",len(RVcatf))
    print("\tCatalogues included:",np.unique(RVcat.CAT.values))
    RVcatf.is_copy=False
    cond=RVcatf.RV==''
    RVcatf=RVcatf.drop(RVcatf.index[cond])
    RVcatf.RV=RVcatf.RV.map(float)
    cond=RVcatf.eRV==''
    RVcatf=RVcatf.drop(RVcatf.index[cond])
    RVcatf.eRV=RVcatf.eRV.map(float)
    print("\tNumber of filtered entries:",len(RVcatf))

    RV=RVcatf
    cats=['hip','tycho2_id']
    for cname in cats:
        RV[cname]=RV[cname].fillna('')
        RV[cname]=RV[cname].map(str)
    dfstr=RV.select_dtypes(['object'])
    RV[dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    RV['hip']=RV['hip'].apply(lambda x:x.replace('.0',''))
    print("\tNumber of RV objects: %d"%len(RV))

    RVCat=pd.DataFrame()
    for col in "hip","tycho2_id":
        sel=RV[RV[col]!=''][[col,"eRV"]].sort_values([col,"eRV"])
        print("\tNumber of total entries for %s: %d"%(col,len(sel)))
        index=sel[col].drop_duplicates().index
        uniq=RV.ix[index]
        print("\tNumber of uniq entries for %s: %d"%(col,len(uniq)))
        RVCat=RVCat.append(uniq)

    
    print("\tTotal number of uniq objects:%d"%len(RVCat))

    dbfile=RV_DIR+"RV.csv"
    print("\tSaving ",dbfile)
    RVCat.to_csv(dbfile,index=False)
    print("\tDone.")

    rv=RVCat

    # =============================================================================
    # MERGING
    # =============================================================================
    print("Merging databases...")
    mdb=database
    cats=['hip','tycho2_id']
    for cname in cats:
        mdb[cname]=mdb[cname].fillna('')
        mdb[cname]=mdb[cname].map(str)
    dfstr=mdb.select_dtypes(['object'])
    mdb[dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
    mdb['hip']=mdb['hip'].apply(lambda x:x.replace('.0',''))

    gaiasrc=database
    cols=["tycho2_id","hip"]
    rvgaia=pd.DataFrame()
    for col in cols[::1],cols[::-1]:
        print("\tMerging by %s..."%col[0])
        result=pd.merge(left=gaiasrc[gaiasrc[col[0]]!=''],
                        right=rv[rv[col[0]]!=''],
                        on=col[0])
        result=result.drop("%s_y"%col[1],1)
        result=result.rename(columns={"%s_x"%col[1]:col[1]})
        print("\tNumber of matchings for %s: %d"%(col[0],len(result)))
        rvgaia=rvgaia.append(result)

    rvgaia=rvgaia.fillna('NULL')
    print("\tNumber of matches: %d"%len(rvgaia))

    dbfile=SRC_DIR+"AstroRV.csv"
    print("\tSaving ",dbfile)
    rvgaia.to_csv(dbfile,index=False)
    print("\tDone.")

if qread:buildCat(gaia,hipparcos,tycho,simbad)
