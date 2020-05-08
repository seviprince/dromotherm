from dromosense import getCsvDatas

step, datas = getCsvDatas("meteo.csv",preview=True)

"""
col 0 : Température air (°C)
col 1 : Température point de rosée (°C)
col 2 : Nature des précipitations
col 3 : Vitesse du vent (m/s)
col 4 : Rayonnement global (W/m2)
col 5 : Rayonnement atmosphérique (W/m2)
"""

"""
Hv Coefficient convectif entre la surface et l'air (fonction affine de la vitesse du vent)
"""
Hv=5.8+4.1*datas[:,3]
B1=Hv*datas[:,0]+datas[:,5]+(1-albedo)*datas[:,4]

