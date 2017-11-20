
# coding: utf-8

# NOTE: Run first the last line

# # Build AstroRV Catalog

# This notebook is aimed to build the AstroRV catalog, a complete compilation of astrometric and radial velocities information of nearby stars. 
# 
# This catalog was originally aimed at studying the interstellar origin probability for objects approaching the Solar System in unbound orbits

# ## Sources

# In[1]:


#GAIA SOURCES
TGAS_DIR="DataGAIA/Data-Gaia/"

#SIMBAD SOURCES
SIMBAD_DIR="DataGAIA/Data-Simbad/"

#TYCHO2/HIPPARCOS SOURCES
HIPTYC_DIR="DataGAIA/Data-Hipparcos/"

#RADIAL VELOCITY SOURCES
RV_DIR="RVGaia/RV/"


# ## Read input catalogs

# ### GAIA

# In[4]:


# =============================================================================
# LECTURA DE ARCHIVOS DE GAIA
# =============================================================================

# Description:

# hip : Hipparcos identifier (int).
# tycho2_id : Tycho 2 identifier (string).
# ref_epoch : Reference epoch (double, Time[Julian Years]), expressed as a Julian Year in TCB.
# ra : Right ascension (double, Angle[deg]). Barycentric right ascension of the source in ICRS.
# ra error : Standard error of right ascension (double, Angle[mas]).
# dec : Declination (double, Angle[deg]). Barycentric declination of the source in ICRS.
# dec error : Standard error of declination (double, Angle[mas]).
# parallax : Parallax (double, Angle[mas]). Absolute barycentric stellar parallax of the source.
# parallax error : Standard error of parallax (double, Angle[mas]).
# pmra : Proper motion in right ascension direction (double, Angular Velocity[mas/year].
# pmra error : Standard error of proper motion in right ascension direction (double, Angular Velocity[mas/year]).
# pmdec : Proper motion in declination direction (double, Angular Velocity[mas/year].
# pmdec error : Standard error of proper motion in declination direction (double, Angular Velocity[mas/year]).
# phot g mean mag : G-band mean magnitude (double, Magnitude[mag]) Mean magnitude in the G band.
# l : Galactic longitude (double, Angle[deg]).
# b : Galactic latitude (double, Angle[deg]).

cols_gaia = ["hip", "tycho2_id", "ra", "ra_error", "dec", "dec_error", "parallax", "parallax_error", "pmra", "pmra_error",         "pmdec", "pmdec_error", "ra_dec_corr", "ra_parallax_corr", "ra_pmra_corr", "ra_pmdec_corr", "dec_parallax_corr",         "dec_pmra_corr", "dec_pmdec_corr", "parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr",         "phot_g_mean_flux", "phot_g_mean_flux_error", "phot_g_mean_mag", "l", "b", "ecl_lon", "ecl_lat"]

for i in range(16):
    filename = TGAS_DIR + "TgasSource_000-000-0" + str(i).zfill(2) + ".csv.gz"
    if i == 0:
        print("Reading", filename)
        gaia = pd.read_csv(filename, usecols=cols_gaia)
    else:
        print("Reading", filename)
        DRx = pd.read_csv(filename, usecols=cols_gaia)
        gaia = gaia.append(DRx)

gaia = pd.DataFrame(gaia)

gaia_hip = gaia[gaia.hip.notnull()]
gaia_tyc = gaia[gaia.tycho2_id.notnull()]

print("Gaia: Subset Hipparcos:", len(gaia_hip))
print("Gaia: Subset Tycho-2:", len(gaia_tyc))
print("Total Gaia:", len(gaia))


# ### HIPPARCOS

# In[5]:


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

# Descartar elementos con paralajes nulos

n1 = len(hipparcos)
#hipparcos.apply(pd.to_numeric, errors="coerce").info()                                # Convierte la información en tipo float

hipparcos["ra_Hip"] = pd.to_numeric(hipparcos["ra_Hip"], errors="coerce")             # Convierte en np.nan los valores no numer.
hipparcos["dec_Hip"] = pd.to_numeric(hipparcos["dec_Hip"], errors="coerce")           # Convierte en np.nan los valores no numer.
hipparcos["parallax_Hip"] = pd.to_numeric(hipparcos["parallax_Hip"], errors="coerce") # Convierte en np.nan los valores no numer.
hipparcos["Vmag_Hip"] = pd.to_numeric(hipparcos["Vmag_Hip"], errors="coerce")         # Convierte en np.nan los valores no numer.
hipparcos.dropna(subset=["ra_Hip", "dec_Hip", "parallax_Hip"], inplace=True)          # Elimina los registros con paralaje nulo
hipparcos["source"]="hipparcos"


n2 = len(hipparcos)
print("Objects discarded:", n1-n2)

len(hipparcos)


# ### TYCHO

# In[6]:


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

# Dividir la cadena de texto en los 3 campos numéricos que componen el ID de Tycho_2 (separados por espacios en blanco)
tycho["a"], tycho["b"], tycho["c"] = tycho["tycho2_id"].str.split().str

# Concatenar los 3 campos numéricos que componen el ID de Tycho_2, separados por guión.
tycho["tycho2_id"] = tycho["a"] + "-" + tycho["b"] + "-" + tycho["c"]

# Borrar campos usados para la conversión
del tycho["a"], tycho["b"], tycho["c"]

# Descartar elementos con paralajes nulos

n1 = len(tycho)
#tycho.apply(pd.to_numeric, errors="coerce").info()                           # Convierte la información en tipo float

tycho["ra_Tyc"] = pd.to_numeric(tycho["ra_Tyc"], errors="coerce")             # Convierte en np.nan los valores no numéricos
tycho["dec_Tyc"] = pd.to_numeric(tycho["dec_Tyc"], errors="coerce")           # Convierte en np.nan los valores no numéricos
tycho["parallax_Tyc"] = pd.to_numeric(tycho["parallax_Tyc"], errors="coerce") # Convierte en np.nan los valores no numéricos
tycho["Vmag_Tyc"] = pd.to_numeric(tycho["Vmag_Tyc"], errors="coerce")         # Convierte en np.nan los valores no numéricos
tycho.dropna(subset=["ra_Tyc", "dec_Tyc", "parallax_Tyc"], inplace=True)      # Elimina los registros con paralaje nulo
tycho["source"]="tycho"

n2 = len(tycho)
print("Objects discarded:", n1-n2)

len(tycho)


# ### SIMBAD

# In[7]:


cols = ["typedident","identifier", "radvel", "coord1(ICRS,J2000/2000)", "plx", "pm", "MagV", "spec.type"]
simbad = pd.read_csv(SIMBAD_DIR+"simbad.csv", usecols=cols, delimiter="|")
simbad = pd.DataFrame(simbad)

# Modificación del ID del catálogo para que quede en formato INTEGER
simbad["hip"] = simbad["typedident"].map(lambda x: str(x)[4:]).astype(float)
del simbad["typedident"]

simbad["plx"] = pd.to_numeric(simbad["plx"], errors="coerce")   

# Cálculo de RA y DEC
# 1. Eliminar espacios en blanco de los lados izquierdo y derecho de cada cadena de texto
simbad["coord1(ICRS,J2000/2000)"] = simbad["coord1(ICRS,J2000/2000)"].str.strip()

# 2. Dividir la cadena de texto en 6 campos (hh:mm:ss para RA y hh:mm:ss para DEC)
simbad["ra_h"], simbad["ra_m"], simbad["ra_s"], simbad["dec_h"], simbad["dec_m"], simbad["dec_s"] =     simbad["coord1(ICRS,J2000/2000)"].str.split(" ").str

# 3. Concatenar los 3 primeros campos mediante la fórmula de conversión de hh:mm:ss a grados para la ascensión recta
simbad["ra_simbad"] = simbad["ra_h"].astype(float)*15 + simbad["ra_m"].astype(float)/60 + simbad["ra_s"].astype(float)/3600

# 4. Concatenar los 3 últimos campos mediante la fórmula de conversión de hh:mm:ss a grados para la declinación
simbad["dec_simbad"] = np.sign(simbad["dec_h"].astype(float)) * (     np.abs(simbad["dec_h"].astype(float)) + simbad["dec_m"].astype(float)/60 + simbad["dec_s"].astype(float)/3600 )

# 5. Borrar campos usados para la conversión
del simbad["coord1(ICRS,J2000/2000)"]
del simbad["ra_h"], simbad["ra_m"], simbad["ra_s"], simbad["dec_h"], simbad["dec_m"], simbad["dec_s"]

# Cálculo del movimiento propio
# 1. Eliminar espacios en blanco de los lados izquierdo y derecho de cada cadena de texto
simbad["pm"] = simbad["pm"].str.strip()

#2. Dividir la cadena de texto en 2 campos (PM_RA y PM_DEC)
simbad["pmra_simbad"], simbad["pmdec_simbad"] = simbad["pm"].str.split(" ").str

# 3. Borrar campos usados para la conversión
del simbad["pm"]

# 4. Elimina los registros con paralaje nulo
n1 = len(simbad)
simbad.dropna(subset=["ra_simbad", "dec_simbad", "plx"], inplace=True)      
n2 = len(simbad)
print("Objects discarded:", n1-n2)

# Formato velocidad radial
simbad["radvel"] = simbad["radvel"].str.strip()                   # Elimina espacios en blanco
simbad["radvel"] = simbad["radvel"].replace("~", np.nan)          # Reemplaza el caracter "~" por Null

# Modificar nombre de algunas columnas
simbad = simbad.rename(columns={"identifier": "name_simbad"})      # Modifica el nombre de la columna de nombres propios
simbad = simbad.rename(columns={"plx": "parallax_simbad"})         # Modifca el nombre de la columna de paralajes
simbad = simbad.rename(columns={"spec.type": "sptype_simbad"})     # Modifca el nombre de la columna de clasificación espectral
simbad = simbad.rename(columns={"radvel": "radial_vel_simbad"})    # Modifica el nombre de la columna de velocidades radiales
simbad = simbad.rename(columns={"MagV": "Vmag_simbad"})            # Modifica el nombre de la columna de magnitud V
simbad["source"]="simbad"

len(simbad)


# ## Merging

# ### Merging with common fields

# In[8]:


# db = Gaia + Hipparcos
gaia["source"] = "gaia"
exclusive_hip = hipparcos[~hipparcos.hip.isin(gaia_hip.hip)]
exclusive_hip.rename(columns=lambda x: x.replace("_Hip", ""), inplace=True)
exclusive_hip["source"] = "hipparcos"
db = pd.concat([gaia, exclusive_hip], axis=0)
print("Aporte Hipparcos:", len(exclusive_hip))

# db = db + Tycho
exclusive_tyc = tycho[~tycho.tycho2_id.isin(gaia_tyc.tycho2_id)]
exclusive_tyc.rename(columns=lambda x: x.replace("_Tyc", ""), inplace=True)
exclusive_tyc["source"] = "tycho"
db = pd.concat([db, exclusive_tyc], axis=0)
print("Aporte Tycho:", len(exclusive_tyc))

# db = db + Simbad
exclusive_simbad = simbad[~simbad.hip.isin(db.hip)]
exclusive_simbad.rename(columns=lambda x: x.replace("_simbad", ""), inplace=True)
exclusive_simbad = exclusive_simbad.loc[:,["hip", "ra", "dec", "parallax", "pmra", "pmdec"]]
exclusive_simbad["source"] = "simbad"
db = pd.concat([db, exclusive_simbad], axis=0)
print("Aporte Simbad:", len(exclusive_simbad))

# db = db + todos los elementos de Simbad [nombre_star, tipo_espectral, velocidad_radial, Vmag]
db = db.merge(simbad.loc[:,["hip", "name_simbad","sptype_simbad","radial_vel_simbad","Vmag_simbad"]],               left_on="hip", right_on="hip", how="left")

# Organizar columnas
cols = ["source"] + cols_gaia + ["Vmag","HenryDraperId"] + ["name_simbad","sptype_simbad","radial_vel_simbad","Vmag_simbad"]
db = db[cols]
db = db.reset_index(drop=True)

# Resultado
print("Tamaño total:", len(db))


# ### Merging with all fields

# In[9]:


n1 = len(gaia)
print("Tamaño inicial de la base de datos alimentada con Gaia:", n1)

database = gaia.merge(hipparcos, left_on="hip", right_on="hip", how="outer")
n2 = len(database)
print("Tamaño de Hipparcos:", len(hipparcos), "Tamaño final de la base de datos maestra:", n2, ". Aporte:", n2-n1)

database = database.merge(tycho, left_on="tycho2_id", right_on="tycho2_id", how="outer")
n3 = len(database)
print("Tamaño de Tycho:", len(tycho), "Tamaño final de la base de datos maestra:", n3, ". Aporte:", n3-n2)

database = database.merge(simbad, left_on="hip", right_on="hip", how="outer")
n4 = len(database)
print("Tamaño de Simbad:", len(simbad), "Tamaño final de la base de datos maestra:", n4, ". Aporte:", n4-n3)

sources=database[['source_x','source_y']].fillna('')
sources.columns=['s1','s2','s3','s4']
database["source"]=sources["s1"].astype(str)+":"+sources["s2"]+":"+sources["s3"]+":"+sources["s4"]
del(database["source_x"])
del(database["source_y"])

database = database.reset_index(drop=True)


# ## Radial velocities

# ### RV Catalogues

# In[15]:


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
print("Building catalogue %s..."%name)
comments=list(range(82))+[83,84]
data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
cats=['HIP']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

#WEB1995
#III/213/catalog
name="Web1995-HIP.csv"
print("Building catalogue %s..."%name)
data[name]=pd.read_csv(srcdir+name)
cats=['HIP']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

#WEB1995-TYC2
#III/213/catalog
name="Web1995-TYC2.csv"
print("Building catalogue %s..."%name)
data[name]=pd.read_csv(srcdir+name)
cats=['TYC2','HIP']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

#GCS2011
#J/A+A/530/A138/catalog
name="GCS2011.tsv"
print("Building catalogue %s..."%name)
comments=list(range(174))+[175,176]
data[name]=pd.read_csv(srcdir+name,sep="|",skiprows=comments)
cats=['HIP']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

#RAVE-DR5
#III/279/rave_dr5
name="RAVE-DR5.tsv"
print("Building catalogue %s..."%name)
comments=list(range(78))+[79,80]
data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
cats=['TYCHO2']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
data[name]["TYC2"]=data[name]["TYCHO2"]
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

#PULKOVO
#III/252/table8
name="Pulkovo.tsv"
print("Building catalogue %s..."%name)
comments=list(range(61))+[62,63]
data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
cats=['HIP']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

#FAMAEY2005
#J/A+A/430/165/tablea1
name="Famaey2005.tsv"
print("Building catalogue %s..."%name)
comments=list(range(118))+[119,120]
data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
cats=['HIP']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

#BB2000
#III/213/catalogue
name="BB2000.csv"
print("Building catalogue %s..."%name)
data[name]=pd.read_csv(srcdir+name)
cats=['TYC2','HIP']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

#MALARODA2012
#III/249/catalog
name="Malaroda2012.csv"
print("Building catalogue %s..."%name)
data[name]=pd.read_csv(srcdir+name)
cats=['TYC2','HIP']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

#GALAH I
#J/MNRAS/465/3203/catal
name="Galah.tsv"
print("Building catalogue %s..."%name)
comments=list(range(54))+[55,56]
data[name]=pd.read_csv(srcdir+name,sep=";",skiprows=comments)
cats=['TYC2']
for cname in cats:
    data[name][cname]=data[name][cname].fillna('')
    data[name][cname]=data[name][cname].map(str)
dfstr=data[name].select_dtypes(['object'])
data[name][dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
print("Number of objects in %s:"%name,len(data[name]))
for cat in cats:
    cond=data[name][cat]!=''
    print("Number of objects in catalogue %s:"%cat,len(data[name][cat][cond]))
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
print("Median error:",med)
print("Number of entries with zero error:",len(df.eRV[df.eRV==0]))
df.eRV[df.eRV==0]=med
#RESULTING SIZE
print("Filtered catalogue:",len(df))
#FILLNA
RVcat=RVcat.append(df.fillna(''))

###################################################
#COMPILING FULL TABLE
###################################################
print("Compiling final catalogue...")
RVcat=RVcat.rename(columns={'TYC2':'tycho2_id','HIP':'hip'})
print("Compiling final catalogue...")
RVcatf=RVcat.reset_index(drop=True)
print("Number of unfiltered entries:",len(RVcatf))
print("Catalogues included:",np.unique(RVcat.CAT.values))
RVcatf.is_copy=False
cond=RVcatf.RV==''
RVcatf=RVcatf.drop(RVcatf.index[cond])
RVcatf.RV=RVcatf.RV.map(float)
cond=RVcatf.eRV==''
RVcatf=RVcatf.drop(RVcatf.index[cond])
RVcatf.eRV=RVcatf.eRV.map(float)
print("Number of filtered entries:",len(RVcatf))

RV=RVcatf
cats=['hip','tycho2_id']
for cname in cats:
    RV[cname]=RV[cname].fillna('')
    RV[cname]=RV[cname].map(str)
dfstr=RV.select_dtypes(['object'])
RV[dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
RV['hip']=RV['hip'].apply(lambda x:x.replace('.0',''))
print("Number of RV objects: %d"%len(RV))

RVCat=pd.DataFrame()
for col in "hip","tycho2_id":
    sel=RV[RV[col]!=''][[col,"eRV"]].sort_values([col,"eRV"])
    print("Number of total entries for %s: %d"%(col,len(sel)))
    index=sel[col].drop_duplicates().index
    uniq=RV.ix[index]
    print("Number of uniq entries for %s: %d"%(col,len(uniq)))
    RVCat=RVCat.append(uniq)

RVCat.to_csv(RV_DIR+"RVCat.csv",index=False)
print("Total number of uniq objects:%d"%len(RVCat))

rv=RVCat


# ## Merging

# In[17]:


#Catalogs
print("GAIA:",len(gaia))
print("HIPPARCOS:",len(hipparcos))
print("TYCHO:",len(tycho))
print("SIMBAD:",len(simbad))
print("RV:",len(rv))


# In[18]:


#Filtering results
database[['source','name_simbad','parallax','parallax_Hip','parallax_Tyc','parallax_simbad','radial_vel_simbad','Vmag_simbad']]


# In[19]:


#Converting hip column into string
mdb=database
cats=['hip','tycho2_id']
for cname in cats:
    mdb[cname]=mdb[cname].fillna('')
    mdb[cname]=mdb[cname].map(str)
dfstr=mdb.select_dtypes(['object'])
mdb[dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
mdb['hip']=mdb['hip'].apply(lambda x:x.replace('.0',''))

#Converting hip column into string
mdb=gaia
cats=['hip','tycho2_id']
for cname in cats:
    mdb[cname]=mdb[cname].fillna('')
    mdb[cname]=mdb[cname].map(str)
dfstr=mdb.select_dtypes(['object'])
mdb[dfstr.columns]=dfstr.apply(lambda x: x.str.strip())
mdb['hip']=mdb['hip'].apply(lambda x:x.replace('.0',''))


# ### Merging GAIA with RV

# In[20]:


#gaiasrc=gaia
gaiasrc=database

cols=["tycho2_id","hip"]
rvgaia=pd.DataFrame()
for col in cols[::1],cols[::-1]:
    print("Merging by %s..."%col[0])
    result=pd.merge(left=gaiasrc[gaiasrc[col[0]]!=''],
                    right=rv[rv[col[0]]!=''],
                    on=col[0])
    result=result.drop("%s_y"%col[1],1)
    result=result.rename(columns={"%s_x"%col[1]:col[1]})
    print("Number of matchings for %s: %d"%(col[0],len(result)))
    rvgaia=rvgaia.append(result)

rvgaia=rvgaia.fillna('NULL')
rvgaia.to_csv(RV_DIR+"RVGaia_all.csv",index=False)
print("Number of matches: %d"%len(rvgaia))


# **Result of merging database Gaia + Hip + Tyc2 + Simbad:**
# 
# Number of matchings for tycho2_id: 230610
# 
# Number of matchings for hip: 37358
# 
# Number of matches: 267968
# 
# **Result of merging database Gaia :**
# 
# Number of matchings for tycho2_id: 210263
# 
# Number of matchings for hip: 25925
# 
# Number of matches: 236188

# In[50]:


",".join(rvgaia.columns)


# In[52]:


rvgaia[['source','CAT','hip','tycho2_id','name_simbad','Vmag_simbad','parallax','parallax_error','parallax_Hip','parallax_error_Hip','radial_vel_simbad','RV','eRV']]


# In[40]:


cond=rvgaia.Vmag_simbad!='NULL'
rvgaia_sel=rvgaia[cond]
rvgaia_sel["Vmag_simbad"]=rvgaia_sel["Vmag_simbad"].map(lambda x:float(x) if x!='~' else 25).astype(float)
#cond=rvgaia.Vmag_simbad<5
#rvgaia_sel=rvgaia_sel[cond]


# In[49]:


rvgaia_sel[rvgaia_sel.Vmag_simbad<2][['name_simbad','parallax','parallax_Hip','parallax_error_Hip','RV','eRV']]


# ## Global Definitions

# In[3]:


import sys
sys.path.append("/home/local-python/lib/python3.5/site-packages")
import numpy as np
import pandas as pd
import collections
pd.options.mode.chained_assignment = None

