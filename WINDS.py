#======================Block 0: Adding Libraries=============================
import sys
import pandas as pd
import numpy as np
from datetime import datetime as dt, timedelta
import matplotlib.pyplot as plt
import pymysql
pymysql.install_as_MySQLdb()
from sqlalchemy import create_engine
#import WINDSfunctionsandclasses13.py

#exec(open("WINDSfunctionsandclasses_July25.py").read())
pathprefix='/home/ecoslacker/Documents/WINDS_Data/'
# sys.path.append(pathprefix)

db=create_engine('mysql://UofABEWINDS:WINDSAWSPort2020@windsdatabase-1.cdzagwevzppe.us-west-1.rds.amazonaws.com:3306/winds_test')

#======================Block 0: Adding Libraries=============================

import WINDSfunctionsandclasses  as wmf

Planting_Array=pd.read_sql('SELECT * from plantings2',con=db)  #Reads all data from mysql db
Soil_Array=pd.read_sql('SELECT * from field_soil_layers',con=db)
Spatial_Array=pd.read_sql('SELECT * from spatial_data',con=db)
WeatherArray=pd.read_sql('SELECT * from weather',con=db)
Average_Weather_Array=pd.read_sql('SELECT * from averageweather',con=db)
Irrigation_Array=pd.read_sql('SELECT * from irrigation_activity',con=db)
ET_dailyArray=pd.read_sql('SELECT * from et_daily',con=db)
Field_Array=pd.read_sql('SELECT * from fields',con=db)
Output_Array=pd.read_sql('SELECT * from output',con=db)
Output_Layer_Array=pd.read_sql('SELECT * from TempLayerOutput',con=db)
ClovisArray=pd.read_sql('SELECT * from Clovis_weather', con=db)
WettingArray=pd.read_sql('SELECT * from Wetted_fractions', con=db)

now = dt.today() 
currentDate = now.strftime("%Y-%m-%d")

#excel_path = pathprefix + 'WINDS guayule_guarNM.xlsx' 
#excel_output_path = pathprefix + 'output_test3.xlsx' 



P = wmf.plantings(Planting_Array) #This sets up a planting object called P with all of the plantings and their associated information
Num_plantings = len(P.DateP) #This is the total number of plantings in the database based on the length of the planting fields.
for i in range(0, Num_plantings): #this loop cycles through all of the plantings
    print(P.RunPlanting[i])    
    if (P.RunPlanting[i] == 1): #checks to see whether a planting is active and should be analyzed
        Run_time = int(P.Max_Days[i])
        P.Set_dates_and_data_sources(i)       
        for kk in range(0, len(Field_Array)):
            field_test = Field_Array.loc[kk]
            print('test and array',field_test['fid'], P.FieldIDP[i])
            if int(field_test['fid']) == int(P.FieldIDP[i]):
                field_index = kk
                print("Made it here")
        for locate in range(1, int(P.NumLocations[i]) + 1): #Loops through the locations in a planting
            print('before soil group letter', i, locate, Num_plantings) #P.PlantingIDP[i]) #, (Spatial_Array['PlantingID']==P.PlantingIDP[i]) & (Spatial_Array['Location']==locate))
            Soil_group_letter = wmf.Find_soil_group(Spatial_Array.loc[(Spatial_Array['PlantingID']==P.PlantingIDP[i]) & (Spatial_Array['Location']==locate)])
            print("i and soil group letter ", i, Soil_group_letter)
            Irrigation_group_letter = wmf.Find_irrigation_group(Spatial_Array.loc[(Spatial_Array['PlantingID']==P.PlantingIDP[i]) & (Spatial_Array['Location']==locate)])
            #P.First day is choosing the date from the last simulation. The last simulation is in output array. Then, it would run from the next day of last simulation to the current date.
            P.First_day = wmf.Check_output_arrays_for_last_DOE_in_previous_output(Output_Array.loc[(Output_Array['pid']==P.PlantingIDP[i]) & (Output_Array['Location']==locate)],
                    Output_Layer_Array.loc[(Output_Layer_Array['pid']==P.PlantingIDP[i]) & (Output_Layer_Array['Location']==locate)], 
                                           P.First_day) #this finds the day that is brought in as day zero from previous simulations
            if P.DateP[i] > dt(2018,1,1) and P.WeatherStationNameP[i] == 'Maricopa':
                Weather2018 = wmf.import_wsdata_online(P.WeatherStationNameP[i], 2018)
                Weather_Array_in = wmf.Create_weather_array(Weather2018.loc[(Weather2018['Date'] > P.DateP[i] + timedelta(int(P.First_day - 1))) & (Weather2018['Date'] <= P.DateP[i] + timedelta(int(P.First_day + P.final_day)))],
                                        Average_Weather_Array.loc[Average_Weather_Array['sid']==P.StationIDP[i]],
                                        P.First_day, P.DateP[i], P.final_day) 
            elif P.WeatherStationNameP[i] == 'Clovis':
                Weather2018 = wmf.Clovis_eto(ClovisArray)
                Weather_Array_in = wmf.Create_weather_array(Weather2018, Average_Weather_Array, P.First_day, P.DateP[i], P.final_day)
            else:
                Weather_Array_in = wmf.Create_weather_array(WeatherArray.loc[(WeatherArray['StationID']==P.StationIDP[i]) & (WeatherArray['Date'] > P.DateP[i] + timedelta(int(P.First_day - 1))) & (WeatherArray['Date'] <= P.DateP[i] + timedelta(int(P.First_day + P.final_day)))],
                                       Average_Weather_Array.loc[Average_Weather_Array['sid']==P.StationIDP[i]],
                                       P.First_day, P.DateP[i], P.final_day) 
            
            Soil_Array_in = wmf.Create_soil_array(Soil_Array.loc[(Soil_Array['FieldID']==P.FieldIDP[i]) & (Soil_Array['SoilGroup']==Soil_group_letter)],
                                    Field_Array.loc[field_index])
            Output_Array_in = Output_Array.loc[(Output_Array['pid']==P.PlantingIDP[i]) & (Output_Array['Location']==locate) & (Output_Array['Days_After_Planting'] == P.First_day)]
            Output_Layer_Array_in = Output_Layer_Array.loc[(Output_Layer_Array['pid']==P.PlantingIDP[i]) & (Output_Layer_Array['Location']==locate) & (Output_Layer_Array['Days_after_planting'] == P.First_day)]                                                                                 
            
            M = wmf.model(Planting_Array.loc[i], Field_Array.loc[field_index],
                  Soil_Array_in, #reads in the soil layer parameters from the layer_data table for the loc ation
                  Weather_Array_in) #reads in the weather values for the season for the station associated with the planting from Weather table

            M.ET_daily(ET_dailyArray.loc[ET_dailyArray['cid']==P.CropID[i]])
            M.SetOutputArrays(P.FirstDate, P.First_day, P.final_day, Output_Array_in)
            M.SetOutputLayersArray(P.First_day, P.final_day, Output_Layer_Array_in)
            M.ParameterizeIrrigation(P.First_day, P.final_day, Irrigation_Array.loc[(Irrigation_Array['pid']==P.PlantingIDP[i]) & (Irrigation_Array['IrrigationGroup']==Irrigation_group_letter)])
            M.Create_output_layer_arrays(int((P.final_day + 1) * (M.Num_layers + 1)), P.final_day, M.Num_layers)

            M.PlantingIDO[0] = P.PlantingIDP[i]
            M.PlantingNameO[0] = P.PlantingName[i]
            M.LocationO[0] = locate
            Name = "Guar_fraction"
            FW_layers = wmf.create_wetted_array(WettingArray, int(P.final_day), M.Num_layers, Name, 500, 1000, 1000, 1000, 1000) 
            
            for j in range(1, int(P.final_day) - 1):
                BL: int = int(M.Bottom_layer[j])
                M.PlantingIDO[j] = P.PlantingIDP[i]
                M.PlantingNameO[j] = P.PlantingName[i]
                M.LocationO[j]  = locate
                for k in range(1, M.Num_layers + 2):
                    M.Fw[k] = M.FW_layers[j][k]
                    M.Fw_y[k]  = M.FW_layers[j - 1][k]
                if int(M.Select_Kcb) == 3:
                    if Kcb_measured > M.Kcb_initial:
                        if Update_Kcb: 
                            M.Kcb[j] = Kcb_measured
                        if Update_Root: 
                            M.Root[j] = Root_measured
                        if Update_Crop_height: 
                            M.Crop_height[j] = Crop_height_measured
                M.Depletion_total[j - 1] =  0
                for k in range(M.Num_layers + 1, BL - 1, -1):
                    M.Percent_depletion[j][k] = (M.FC[k] - M.WC[j-1][k]) /(M.FC[k] - M.PWP[k]) * 100
                    M.Depletion_total[j - 1] = M.Depletion_total[j - 1] + (M.FC[k] - M.WC[j - 1][k]) * M.dz[k] *  M.Fw_y[k]
                M.Percent_depletion_total[j-1] = M.Depletion_total[j-1] / M.Total_TAW_Fw[j][BL] * 100      
                if M.Neglect_upper_layer_in_depletion_calc == 1:
                    M.Percent_depletion_total[j - 1] = (M.Depletion_total[j - 1] - (M.FC[M.Num_layers + 1] - M.WC[j - 1][M.Num_layers + 1]) * M.dz[M.Num_layers + 1] * M.Fw_y[M.Num_layers + 1]) / (M.Total_TAW_Fw[j][BL] - M.TAW[M.Num_layers + 1] * M.Fw_y[M.Num_layers + 1]) * 100
                if M.adjust_P == 1:
                    if M.Ksat[M.Num_layers + 1] * 100 / 24 > 50:
                        p_adj = M.pTable22 + 0.1
                    elif (M.Ksat[M.Num_layers + 1] * 100 / 24) > 25:
                        p_adj = M.pTable22 +0.05
                    elif Ksat[M.Num_layers + 1] * 100 / 24 < 0.5:
                        p_adj = M.pTable22 -  0.1
                    elif Ksat[M.Num_layers + 1] * 100 / 24 < 2.5:
                        p_adj = M.pTable22 - 0.05
                    else:
                        p_adj = M.pTable22
                    if j == 1: 
                        M.ET_trans[j - 1] = 5
                    M.Ps[j] = p_adj + 0.04 * (5 - (M.ET_trans[j - 1]))
                else:
                    M.Ps[j] = M.pTable22
                
                if int(M.Select_Kcb) == 6:
                    M.Ps[j] = 0.75
                
                #this is for showing the dynamic simulation 
                
                M.RAWsum[j] = M.TAWsum[j] * M.Ps[j]
                M.Ks_water_total[j] = (M.Total_TAW_Fw[j][BL] - M.Depletion_total[j - 1]) / ((1 - M.Ps[j]) * M.Total_TAW_Fw[j][BL])
                if M.Ks_water_total[j] < 0: 
                    M.Ks_water_total[j] = 0
                M.Ks[j] = M.Ks_water_total[j]
                        
                if M.Ks[j] > 1:
                    M.Ks[j] = 1

                if M.Salinity_simulation == 1:
                    M.Ks_salt[j] = 1 - M.b_sal / (100) * (M.ECe_ave_effective[j - 1] - M.ECe_thresh)
                    if M.Ks_salt[j] > 1: 
                        M.Ks_salt[j] = 1
                    if M.Ks_salt[j] < 0: 
                        M.Ks_salt[j] = 0
                else:
                    M.Ks_salt[j] = 1
        
                if M.TEW > 0:
                   
                    if M.RH_min[j] < 20:
                        RH_used = 20
                    else:
                        RH_used = M.RH_min[j]

                    M.Kcmax[j] = 1.2 + (0.04 * (M.WS[j - 1] - 2) - 0.004 * (RH_used - 45)) * (M.Crop_height[j] / 3) ** 0.3
                    
                    if M.Kcb[j] + 0.05 > M.Kcmax[j]:
                        M.Kcmax[j] = M.Kcb[j] + 0.05
                    
                    if M.Kcb[j] < 0.15:
                        one_minus_c = 1
                    else:
                        one_minus_c = 1 - ((M.Kcb[j] - 0.15) / (M.Kcmax[j] - 0.15)) ** (1 + M.Crop_height[j] * 0.5)
                    if one_minus_c < 0.01:
                        one_minus_c = 0.01
                    if M.Kcb[j] <= 0.15:
                        one_minus_c = 1
                    else:
                        one_minus_c = 1 - ((M.Kcb[j] - 0.15) / (M.Kcmax[j] - 0.15)) ** (1 + M.Crop_height[j] * 0.5)
                   
                    
                    if int(M.Select_Kcb) == 6:
                        if j >= M.dev:
                            one_minus_c = M.ETone_minus_c_daily[j]
                        else:
                            one_minus_c = 1
                    Active_depletion = (M.FC[M.Num_layers + 1] - M.WC[j - 1][M.Num_layers + 1]) * M.dz[M.Num_layers + 1] * 1000
                    if Active_depletion > M.REW:
                        M.Kr[j] = (M.TEW - Active_depletion) / (M.TEW - M.REW)
                        if M.Kr[j] > 1:
                            M.Kr[j] = 1
                        elif M.Kr[j] < 0:
                            M.Kr[j] = 0
                    else:
                        M.Kr[j] = 1
                    
                    if (M.Kcmax[j] * one_minus_c) < ((M.Kcmax[j] - M.Kcb[j]) * M.Kr[j] * M.Fw[M.Num_layers + 1]):
                        M.Ke[j] = M.Kcmax[j] * one_minus_c
                    else:
                        M.Ke[j] = (M.Kcmax[j] - M.Kcb[j]) * M.Kr[j] * M.Fw[M.Num_layers + 1]
                
                           
                else:
                    M.Ke[j] = 0

                if M.No_stress_reduction == 1:
                    M.Ks[j] = 1#
                    M.Ks_salt[j] = 1#
                
                M.Ks_salt[j] = 1
                M.ET_trans[j] = M.Kcb[j] * M.Ks[j] * M.Ks_salt[j] * M.ET_0[j]
                M.ET_pot[j] = M.Kcb[j] * M.ET_0[j]
                M.E[j] = M.Ke[j] * M.ET_0[j]
                
                   
                if M.No_ET == 1:
                    M.ET_trans[j] = 0
                    M.E[j] = 0
                
                M.ET_actual_cum[j] = M.ET_actual_cum[j - 1] + M.ET_trans[j]
                M.ETcum[j] = M.ETcum[j - 1] + M.ET_pot[j]
                M.ET[j] = M.ET_trans[j] + M.E[j]    
                M.Depletion_total[j] = M.Depletion_total[j - 1] + (M.ET_trans[j] + M.E[j]) / 1000
                M.Percent_depletion_total[j] = M.Depletion_total[j] / M.Total_TAW_Fw[j][BL] * 100
                M.AW[j] = M.TAWsum[j] - M.Depletion_total[j] - M.RAWsum[j]  
                if M.Rainfall_simulation == 1:
                    if M.rainfall[j] > 0:
                        if M.Rain_infiltration_calculation == 1:
                            M.Rain_Infilt[j] = 1000 * M.wetting_front(j)
                        else:
                            M.Rain_Infilt[j] = M.rainfall[j] * M.Rain_partition
                        Rain_happened = 0
                    else:
                        M.Rain_Infilt[j] = 0
                
                if M.No_infiltration == 1:
                    M.Irrigation[j] = 0
                    M.Rain_Infilt[j] = 0
                else:
                    M.Irrigation[j] = M.Irrigation_depth[j]
#         Allocate transpiration between layers based on depth of roots and actual depletion
                sum_ET = 0
                for k in range(

                        BL, M.Num_layers + 2):
                    sum_ET = sum_ET + M.ETfrac[k]
        
                if M.No_ET_frac_Adjustment == 1:
                    for k in range(1, M.Num_layers + 2):
                        M.Act_frac[j][k] = M.ETfrac[k] / sum_ET
                else:
                    Remaining = 1
                    Sum_act_frac = 0
                    
                    for k in range(M.Num_layers + 1, BL - 1, -1):
                        M.RAW[k] = M.Ps[j] * M.TAW[k]
                        if M.Fw[k] < M.Fw_y[k]:
                            M.Depletion[j - 1][k] = M.Depletion[j - 1][k] * M.Fw[k] / M.Fw_y[k]
                        if (M.Depletion[j - 1][k] > M.RAW[k] * M.Fw[k]) and (M.Depletion[j - 1][k] < M.TAW[k] * M.Fw[k]):
                            M.Act_frac[j][k] = M.ETfrac[k] * (M.TAW[k] * M.Fw[k] - M.Depletion[j - 1][k]) / ((M.TAW[k] - M.RAW[k]) * M.Fw[k])/sum_ET
                        elif M.Depletion[j - 1][k] > M.TAW[k] * M.Fw[k]:
                            M.Act_frac[j][k] = 0
                        else:
                            M.Act_frac[j][k] = M.ETfrac[k]/sum_ET
                        Sum_act_frac = Sum_act_frac + M.Act_frac[j][k]
                    
                    if Sum_act_frac > 0:
                        for k in range(M.Num_layers + 1, BL - 1, -1):
                            M.Act_frac[j][k] = M.Act_frac[j][k] / Sum_act_frac

                te = np.zeros(M.Num_layers + 2)
                kh = np.zeros(M.Num_layers + 2)
                Keff = np.zeros(M.Num_layers + 2)
                
                M.E_wet[j] = M.E[j] * M.Fw[M.Num_layers + 1]
                
                   
                M.wet_Rain_Infilt[j] = M.Rain_Infilt[j] * M.Fw[M.Num_layers + 1]
                if M.Water_table_simulation == 1:
                    M.Equilibrium_Max[j] = M.Equilibrium_Max[j-1]
                else:
                    M.Equilibrium_Max[j] = 0
                    
                Equil_Max = int(M.Equilibrium_Max[j])
                M.TDW[j] = 0
                
                Depthcalc = 0
                
                for k in range(M.Num_layers + 1, Equil_Max, -1): # Equil_Max is the upper layer number of the water table. 
                   
                    M.Active_WC[k] = M.WC[j - 1][k] # M.Active is the water content ###Probably not needed for irrigation. If not, replace it with M.WC
                    
                    if k == M.Num_layers + 1:
                        M.WC[j][k] = M.Active_WC[k] + (M.Irrigation[j] - M.E_wet[j] - M.ET_trans[j] * M.Act_frac[j][k] + M.wet_Rain_Infilt[j]) / M.Fw[k]  / 1000 / M.dz[k]
                    else:
                        M.WC[j][k] = M.Active_WC[k] + (M.Infilt[j][k + 1] - M.ET_trans[j] / 1000 * M.Act_frac[j][k]) / M.dz[k] / M.Fw[k]
                    
                    if M.WC[j][k] > M.FC[k] and M.Turn_off_field_capacity_restriction == 0:
                        M.Infilt[j][k] = (M.WC[j][k] - M.FC[k]) * M.dz[k] * M.Fw[k]
                        M.WC[j][k] = M.FC[k]
                    elif M.WC[j - 1][k - 1] >= M.ts[k - 1]:
                        M.WC[j - 1][k - 1] = M.ts[k - 1]
                        M.Infilt[j][k] = 0
                    elif M.Active_WC[k] > M.PWP[k] and k != 1:
                        te[k] = wmf.Calculate_te(M.Active_WC[k], M.tr[k], M.ts[k])
                        te[k - 1] = wmf.Calculate_te(M.WC[j - 1][k - 1], M.tr[k - 1], M.ts[k - 1])
                        kh[k] = wmf.Calculate_k(te[k], M.mv[k], M.Ko[k], M.Lv[k])
                        kh[k - 1] = wmf.Calculate_k(te[k - 1], M.mv[k - 1], M.Ko[k - 1], M.Lv[k - 1])
                        Keff1 = (M.dz[k] + M.dz[k - 1]) / (M.dz[k] / kh[k] + M.dz[k - 1] / kh[k - 1])
                        Keff2 = (kh[k] + kh[k - 1]) / 2
                        if te[k] < te[k - 1] or M.WC[j][k] < M.FC[k]:
                            Keff[k] = Keff1
                        else:
                            Keff[k] = Keff2
                        
                        M.hc[j][k] = -((te[k] ** (-1 / M.mv[k]) - 1) ** (1 / M.nv[k])) / M.av[k] / 100
                        M.hc[j][k - 1] = -((te[k - 1] ** (-1 / M.mv[k - 1]) - 1) ** (1 / M.nv[k - 1])) / M.av[k - 1] / 100
                        if k == 1 and Seal_bottom == 1:
                            M.Infilt[j][k] = 0
                        else:
                            M.Infilt[j][k] = (M.hc[j][k] - M.hc[j][k - 1] + (M.Zave[k] - M.Zave[k - 1])) / (M.Zave[k] - M.Zave[k - 1]) * Keff[k] * M.Fw[k]
                        
                        M.WC[j][k] = M.WC[j][k] - M.Infilt[j][k] / M.dz[k] / M.Fw[k]
                    else:
                        M.Infilt[j][k] = 0
                  
                
                    M.Depletion[j][k] = (M.FC[k] - M.WC[j][k]) * M.dz[k] * M.Fw[k]
                    M.Percent_depletion[j][k] = (M.FC[k] - M.WC[j][k]) / (M.FC[k] - M.PWP[k]) * 100
                    if k == 1:
                        M.WC[j][0] = M.WC[j][1]
                    M.TDW[j] = M.TDW[j] + M.WC[j][k]*M.dz[k]
                    Depthcalc = Depthcalc + M.dz[k]
                M.VWCave[j] = M.TDW[j] / Depthcalc
                        
                if M.Neglect_upper_layer_in_depletion_calc == 1:
                    M.Percent_depletion_total[j] = (M.Depletion_total[j] - (M.FC[M.Num_layers + 1] - M.WC[j][M.Num_layers + 1]) * M.dz[M.Num_layers + 1] * M.FW_layers[j - 1][M.Num_layers + 1]) / (M.Total_TAW_Fw[j][BL] - M.TAW[M.Num_layers + 1] * M.Fw[M.Num_layers + 1]) * 100

                if M.Salinity_simulation == 1:
                    Check = 1
                    for k in range(M.Num_layers + 1, 0, -1):
                        M.Active_EC[k] = M.EC[j - 1][k]

                    for k in range(M.Num_layers + 1, 0, -1):
                     
                    #Calculate salinity addition by wastewater or manure
                    
                        if k == M.Num_layers + 1:
                            mwaste = 0.5 * M.Waste_sal[j] / 10 / 640
                        elif k == M.Num_layers:
                            mwaste = 0.5 * M.Waste_sal[j] / 10 / 640
                        else:
                            mwaste = 0
                        
                        if k == M.Num_layers + 1:
                            if (M.Irrigation[j] + M.wet_Rain_Infilt[j]) / 1000 / M.Fw[k] > (M.FC[k] - M.Active_WC[k]) * M.dz[k]:
                                M.EC[j][k] = (M.Irrigation[j] / 1000 / M.Fw[k] * M.Irr_Sal[j] + mwaste + M.Active_EC[k] * M.dz[k] * M.Active_WC[k]) / ((M.Irrigation[j] + M.wet_Rain_Infilt[j]) / 1000 / M.Fw[k] + M.dz[k] * M.Active_WC[k] - (M.E_wet[j] + M.ET_trans[j] * M.Act_frac[j][k]) / M.Fw[k] / 1000)
                                M.Active_EC[k] = M.EC[j][k]
                                Check = 1
                            elif M.Infilt[j][k] < 0:
                                M.EC[j][k] = (-M.Infilt[j][k] / M.Fw[k] * M.Active_EC[M.Num_layers] + mwaste + M.Active_EC[k] * M.dz[k] * M.Active_WC[k]) / (M.dz[k] * M.WC[j][k])
                                Check = 0
                            else:
                                M.EC[j][k] = (M.Irrigation[j] / 1000 / M.Fw[k] * M.Irr_Sal[j] + mwaste - M.Infilt[j][k] * M.Active_EC[k] / M.Fw[k] + M.Active_EC[k] * M.dz[k] * M.Active_WC[k]) / (M.dz[k] * M.WC[j][k])
                                Check = 0
                            
                            M.ECe[j][k] = M.EC[j][k] * M.WC[j][k] / M.ts[k]
                        else:
                            if Check == 1 and M.Infilt[j][k] > 0:
                                M.EC[j][k] = (M.Infilt[j - 1][k] / M.Fw[k] * M.Active_EC[k + 1] + mwaste + M.Active_EC[k] * M.dz[k] * M.Active_WC[k]) / (M.Infilt[j - 1][k] / M.Fw[k] + M.dz[k] * M.Active_WC[k] - M.ET_trans_wet[j] * M.Act_frac[j][k] / M.Fw[k] / 1000)
                                M.Active_EC[k] = M.EC[j][k]
                                Check = 1
                                if k == 1:
                                    M.EC_leach_eqn[j] = 1
                            elif M.Infilt[j - 1][k] <= 0 and M.Infilt[j][k] > 0:
                                M.EC[j][k] = (M.Infilt[j - 1][k] * M.Active_EC[k] / M.Fw[k] + mwaste - M.Infilt[j][k] * M.Active_EC[k] / M.Fw[k] + M.Active_EC[k] * M.dz[k] * M.WC[j - 1][k]) / (M.dz[k] * M.WC[j][k])
                            elif M.Infilt[j - 1][k] <= 0 and M.Infilt[j][k] <= 0:
                                M.EC[j][k] = (M.Infilt[j - 1][k] * M.Active_EC[k] / M.Fw[k] + mwaste - M.Infilt[j][k] * M.Active_EC[k - 1] / M.Fw[k] + M.Active_EC[k] * M.dz[k] * M.WC[j - 1][k]) / (M.dz[k] * M.WC[j][k])
                            elif M.Infilt[j - 1][k] > 0 and M.Infilt[j][k] <= 0:
                                M.EC[j][k] = (M.Infilt[j - 1][k] * M.Active_EC[k + 1] / M.Fw[k] + mwaste - M.Infilt[j][k] * M.Active_EC[k - 1] / M.Fw[k] + M.Active_EC[k] * M.dz[k] * M.WC[j - 1][k]) / (M.dz[k] * M.WC[j][k])
                            else:
                                M.EC[j][k] = (M.Infilt[j - 1][k] * M.Active_EC[k + 1] / M.Fw[k] + mwaste - M.Infilt[j][k] * M.Active_EC[k] / M.Fw[k] + M.Active_EC[k] * M.dz[k] * M.WC[j - 1][k]) / (M.dz[k] * M.WC[j][k])
            
                            if M.EC[j][k] > M.Max_S:
                                M.EC[j][k] = M.Max_S
                            M.ECe[j][k] = M.EC[j][k] * M.WC[j][k] / M.ts[k]    
                            Act_frac_sum = 0
                            M.Total_salt[j] = 0
                            for k in range(M.Num_layers, BL - 1, -1):
                                M.ECe_ave_effective[j] = M.ECe_ave_effective[j] + M.Act_frac[j][k] * M.ECe[j][k] * M.FC[k] / M.WC[j][k]
                                Act_frac_sum = Act_frac_sum + M.Act_frac[j][k]
                            if Act_frac_sum == 0: 
                                Act_frac_sum = 1
                            M.ECe_ave_effective[j] = M.ECe_ave_effective[j] / Act_frac_sum
                            for k in range(1, M.Num_layers + 2):
                                M.Total_salt[j] = M.Total_salt[j] + M.EC[j][k] * 640 * M.WC[j][k] * M.dz[k] * M.Fw[k]
                            M.EC[j][0] = M.EC[j][1]
                        
                if M.Nitrogen_simulation == 1:
                    
                    Weighted_N = np.zeros(M.Num_layers + 2)
                    
                    M.Total_nit[j - 1] = 0
                    for k in range(1, M.Num_layers + 2):
                        M.Total_nit[j - 1] = M.Total_nit[j - 1] + M.N[j - 1][k] * M.WC[j - 1][k] * M.dz[k] * M.Fw[k]
                    for k in range(1, M.Num_layers + 2):
                        M.Active_n[k] = (M.WC[j - 1][k] * M.N[j - 1][k]) / M.WC[j - 1][k]
                
                    M.Total_Min[j] = 0
                    M.Total_Den[j] = 0
                    M.Total_Fer[j] = 0
                    M.Total_Upt[j] = 0
                    M.Total_nit_accum[j] = 0
                    M.Total_nit[j] = 0
            
                    Total_N = 0
                    
                    N_total_weight = 0
                    Total_weighting = 0
                    Frac_total_weight = 0
                    Total_frac_weighting = 0
                    
                    for k in range(BL, M.Num_layers + 2):
                        if M.WC[j][k] > M.PWP[k]:
                            Weighted_N[k] = M.N_soil[j - 1][k] * M.Act_frac[j][k]
                            N_total_weight = Weighted_N[k] + N_total_weight
                            Total_weighting = Total_weighting + M.Act_frac[j][k]
                        else:
                            Weighted_N[k] = 0
                    
                    if Total_weighting == 0:
                        Total_weighting = 1
                    M.N_ave[j] = N_total_weight / Total_weighting
                    
                    average_Fw = (M.Fw[M.Num_layers] + M.Fw[M.Num_layers + 1]) / 2
                    
                    mupt_optimal = M.N_upt[j] / 10 * M.Frac_NO3 * average_Fw #'Adjusts uptake for units and fraction taken up as nitrate from soil
                    
                    M.N_max[j] = mupt_optimal * (1 + M.Frac_greater_saturation)
                    
            #'Calculate Km for Michaelis-Menton, this could be specified in worksheet
            
                    if M.Constant_Km == 1:
                        Km = M.Km_daily[j]
                    elif M.Constant_Km == 2:
                        if M.N_upt[j] > 0:
                            Km = (M.N_soil_optimal - M.Nmin) * (M.N_max[j] / mupt_optimal - 1)
                        else:
                            Km = 100000
                    else:
                        Km = M.Km
                    
                    M.mupt_entire_profile[j] = M.N_max[j] * (M.N_ave[j] - M.Nmin) / (Km + M.N_ave[j] - M.Nmin)
                    
                    M.N_optimal_low[j] = mupt_optimal * (1 - M.Range_frac / 2)
                    M.N_optimal_high[j] = mupt_optimal * (1 + M.Range_frac / 2)
                    
                    if M.mupt_entire_profile[j] > 0:
                        if M.mupt_entire_profile[j] < M.N_optimal_low[j]:
                            M.Ks_nit[j] = 1 - M.KN_low[j] * (M.N_optimal_low[j] - M.mupt_entire_profile[j]) / M.N_optimal_low[j]
                        elif M.mupt_entire_profile[j] < M.N_optimal_high[j]:
                            M.Ks_nit[j] = 1
                        else:
                            M.Ks_nit[j] = 1 - M.KN_high[j] * (M.mupt_entire_profile[j] - M.N_optimal_high[j]) / M.N_optimal_high[j]
                    else:
                        M.Ks_nit[j] = 1#
                    
                    if M.Ks_nit[j] < 0:
                        M.Ks_nit[j] = 0
                     
                    Weighted_N_ET = 0
                    
                    for k in range(BL, M.Num_layers + 2):
                        if M.WC[j][k] > M.PWP[k]:
                            Weighted_N_ET = Weighted_N_ET + M.ETfrac[k] * M.N_soil[j - 1][k]
                    
                    for k in range(M.Num_layers + 1, 0, -1):
            
                        if (k > BL - 1) and (M.WC[j][k] > M.PWP[k]):
                            mupt = M.mupt_entire_profile[j] * M.ETfrac[k] * M.N_soil[j - 1][k] / Weighted_N_ET
                        else:
                            mupt = 0
                        M.Uptake[j][k] = mupt
            
            #'Calculate mineralization rate
            
                        if M.WC[j - 1][k] < M.PWP[k]:
                            fmnq = 0
                        elif M.WC[j - 1][k] < M.qlow:
                            fmnq = (M.WC[j - 1][k] - M.PWP[k]) / (M.qlow - M.PWP[k])
                        elif M.WC[j - 1][k] < M.qhigh:
                            fmnq = 1
                        else:
                            fmnq = 0.6 + 0.4 * (M.ts[k] - M.WC[j - 1][k]) / (M.ts[k] - M.qhigh)
                        
                        Gmn = M.Kmnl * M.fmntemp[j][k] * fmnq * M.bd[k] * M.Org[k]
                        
                        mMin = Gmn * M.dz[k] * M.FW_layers[j - 1][k]
                        
            #'Calculate denitrification rate
            
                        if M.WC[j - 1][k] < 0.6 * M.ts[k]:
                            fdenq = 0
                        else:
                            fdenq = ((M.WC[j - 1][k] - 0.6 * M.ts[k]) / (M.ts[k] - 0.6 * M.ts[k])) ** 2
                        
                        Gden = M.Kden * M.fmntemp[j][k] * fdenq * M.N[j - 1][k] * M.WC[j - 1][k] * np.exp(-M.alpha * M.Z_top_ave[k] * 100)
                        
                        mDen = Gden * M.dz[k] * M.FW_layers[j - 1][k]
                        
            #'Calculate fertilizer addition
            
                        if k == M.Num_layers + 1:
                            if M.Dont_fertilize_evap_layer == 1:
                                mfer = 0
                            else:
                                mfer = 0.5 * M.Fert[j] / 10
                        elif k == M.Num_layers:
                            if M.Dont_fertilize_evap_layer == 1:
                                mfer = M.Fert[j] / 10
                            else:
                                mfer = 0.5 * M.Fert[j] / 10
                        else:
                            mfer = 0
                        
            #'Make calculation of nitrogen concentrations
                        
                        if k == M.Num_layers + 1:
                                if (M.Irrigation[j] + M.wet_Rain_Infilt[j]) / 1000 / M.Fw[k] > (M.FC[k] - M.Active_WC[k]) * M.dz[k]:
                                    M.N[j][k] = (M.Irrigation[j] / 1000 / M.Fw[k] * M.Irr_Nit[j] + mMin - mDen + mfer - mupt / M.Fw[k] + M.Active_n[k] * M.dz[k]* M.Active_WC[k]) / ((M.Irrigation[j] + M.wet_Rain_Infilt[j]) / 1000 / M.Fw[k] + M.dz[k] * M.Active_WC[k] - (M.E_wet[j] + M.ET_trans_wet[j] * M.Act_frac[j][k]) / M.Fw[k] / 1000)
                                    M.Active_n[k] = M.N[j][k]
                                    Check = 1
                                elif M.Infilt[j][k] < 0:
                                    M.N[j][k] = (-M.Infilt[j][k] / M.Fw[k] * M.Active_n[M.Num_layers] + mMin - mDen + mfer - mupt / M.Fw[k] + M.Active_n[k] * M.dz[k] * M.Active_WC[k]) / (M.dz[k] * M.WC[j][k])
                                    Check = 0
                                else:
                                    M.N[j][k] = (M.Irrigation[j] / 1000 / M.Fw[k] * M.Irr_Nit[j] + mMin - mDen + mfer - mupt / M.Fw[k] - M.Infilt[j][k] * M.Active_n[k] / M.Fw[k] + M.Active_n[k] * M.dz[k] * M.Active_WC[k]) / (M.dz[k] * M.WC[j][k])
                                    Check = 0
                                
                        else:
                                if Check == 1 and M.Infilt[j][k] > 0:
                                    M.N[j][k] = (M.Infilt[j][k + 1] / M.Fw[k] * M.Active_n[k + 1] + mMin - mDen + mfer - mupt / M.Fw[k] + M.Active_n[k] * M.dz[k] * M.Active_WC[k]) / (M.Infilt[j][k + 1] / M.Fw[k] + M.dz[k] * M.Active_WC[k] - M.ET_trans_wet[j] * M.Act_frac[j][k] / M.Fw[k] / 1000)
                                    Check = 1
                                    M.Active_n[k] = M.N[j][k]
                                    if k == 1:
                                        M.N_leach_eqn[j] = 1
                                elif M.Infilt[j][k + 1] <= 0 and M.Infilt[j][k] > 0:
                                    M.N[j][k] = (M.Infilt[j][k + 1] * M.Active_n[k] / M.Fw[k] + mMin - mDen + mfer - mupt / M.Fw[k] - M.Infilt[j][k] * M.Active_n[k] / M.Fw[k] + M.Active_n[k] * M.dz[k] * M.WC[j - 1][k]) / (M.dz[k] * M.WC[j][k])
                                elif M.Infilt[j][k + 1] <= 0 and M.Infilt[j][k] <= 0:
                                    M.N[j][k] = (M.Infilt[j][k + 1] * M.Active_n[k] / M.Fw[k] + mMin - mDen + mfer - mupt / M.Fw[k] - M.Infilt[j][k] * M.Active_n[k - 1] / M.Fw[k] + M.Active_n[k] * M.dz[k] * M.WC[j - 1][k]) / (M.dz[k] * M.WC[j][k])
                                elif M.Infilt[j][k + 1] > 0 and M.Infilt[j][k] <= 0:
                                    M.N[j][k] = (M.Infilt[j][k + 1] * M.Active_n[k + 1] / M.Fw[k] + mMin - mDen + mfer - mupt / M.Fw[k] - M.Infilt[j][k] * M.Active_n[k - 1] / M.Fw[k] + M.Active_n[k] * M.dz[k] * M.WC[j - 1][k]) / (M.dz[k] * M.WC[j][k])
                                else:
                                    M.N[j][k] = (M.Infilt[j][k + 1] * M.Active_n[k + 1] / M.Fw[k] + mMin - mDen + mfer - mupt / M.Fw[k] - M.Infilt[j][k] * M.Active_n[k] / M.Fw[k] + M.Active_n[k] * M.dz[k] * M.WC[j - 1][k]) / (M.dz[k] * M.WC[j][k])
                        M.N_soil[j][k] = M.N[j][k] * M.WC[j][k] / ((1 - M.ts[k]) * 2.65)
                        M.kgha[j][k] = M.N_soil[j][k] * (1 - M.ts[k]) * 2.65 * M.dz[k] * 10
                        M.Total_Min[j] = mMin * M.Fw[k] + M.Total_Min[j]
                        M.Total_Den[j] = mDen * M.Fw[k] + M.Total_Den[j]
                        M.Total_Fer[j] = mfer * M.Fw[k] + M.Total_Fer[j]
                        M.Total_Upt[j] = mupt * M.Fw[k] + M.Total_Upt[j]
                        M.Net_accumulation[j][k] = (mMin - mDen + mfer - mupt)

                    M.Cum_Min[j] = M.Cum_Min[j - 1] + M.Total_Min[j] * 10
                    M.Cum_Den[j] = M.Cum_Den[j - 1] - M.Total_Den[j] * 10
                    M.Cum_Fer[j] = M.Cum_Fer[j - 1] + M.Total_Fer[j] * 10
                    M.Cum_Upt[j] = M.Cum_Upt[j - 1] - M.Total_Upt[j] * 10
                    M.Cum_Drain[j] = M.Cum_Drain[j - 1] - M.Infilt[j][1] * M.N[j][1] * 10
                    M.Cum_Irr[j] = M.Cum_Irr[j - 1] + M.Irrigation[j] / 1000 * M.Irr_Nit[j] * 10
                    M.Cum_Total[j] = M.Cum_Min[j] + M.Cum_Den[j] + M.Cum_Fer[j] + M.Cum_Upt[j] + M.Cum_Drain[j] + M.Cum_Irr[j]
                    for k in range(1, M.Num_layers + 2):
                        M.Total_nit_accum[j] = M.Total_nit_accum[j] + M.Net_accumulation[j][k]
                        M.Total_nit[j] = M.Total_nit[j] + M.N[j][k] * M.WC[j][k] * M.dz[k] * M.Fw[k]
                    M.N[j][0] = M.N[j][1]
                    
            temp_output = pd.DataFrame(columns = ['Date','Days_after_planting', 'DOY', 'pid', 'plantname', 'Location', 'NDVI', 'NDVI_regression', 'Kcb_NDVI', 'Kcb_calculated', 'Kcb_used',
                                                  'one_minus_c', 'Crop_height', 'Root_depth', 'Ke', 'Evap', 'Transpiration', 'ET', 'Pottranspiration', 'CumulativePotentialET',
                                                  'CumulativeactualET', 'ReferenceET', 'Ky', 'Depletion_total', 'Percent_depletion_total', 'Irrigation', 'P', 'Rainfall', 'Rain_infilt',
                                                  'Field_capacity_ave', 'Permanent_wilting_point_ave', 'TAWsum', 'RAWsum', 'Available_water', 'Available_water_depth_ave', 'Volumetric_water_content_ave', 'Total_depth_water',
                                                  'Watertableelevation', 'HeightEqu', 'NumberEqu', 'WaterstressKs', 'IrrigationSaltAdded', 'WasteSaltAdded', 'SeepageSaltLost', 'SumOfFluxes', 
                                                  'TotalSalt', 'DifferenceInSalt', 'SaltStress', 'NitStress', 'TotalIrrN', 'TotalDrainN', 'SumReactionsN', 'SumFluxReactionsN', 'TotalMin',
                                                  'TotalDen', 'TotalFer', 'TotalUpt', 'UptakeReq', 'TotalNit', 'TotalNitAccum', 'NitAve', 'NOptimalLow', 'NOptimalHigh', 'CumMin', 'CumDen', 
                                                  'CumFer', 'CumUpt', 'CumIrrN', 'CumDrnN', 'CumChangeN'])                 

            for iii in range(1, P.final_day - 1):
                M.datelist[iii] = dt.fromordinal(M.datelist[iii].toordinal())

            temp_output['Date']=M.datelist[1:P.final_day - 1]
            temp_output['Days_after_planting']=M.DOE[1:P.final_day - 1]
            temp_output['DOY']=M.DOY[1:P.final_day - 1]
            temp_output['pid']=M.PlantingIDO[1:P.final_day - 1]
            temp_output['plantname']=M.PlantingNameO[1:P.final_day - 1]
            temp_output['Location']=M.LocationO[1:P.final_day - 1]
            temp_output['NDVI'] = M.NDVI[1:P.final_day - 1]
            temp_output['NDVI_regression'] = M.NDVI_regression[1:P.final_day - 1]
            temp_output['Kcb_NDVI'] = M.Kcb_NDVI[1:P.final_day - 1]
            temp_output['Kcb_calculated'] = M.Kcb_calculated[1:P.final_day - 1]
            temp_output['Kcb_used'] = M.Kcb[1:P.final_day - 1]
            temp_output['one_minus_c'] = M.one_minus_c[1:P.final_day - 1]
            temp_output['Crop_height'] = M.Crop_height[1:P.final_day - 1]
            temp_output['Root_depth'] = M.Root[1:P.final_day - 1]
            temp_output['Ke'] = M.Ke[1:P.final_day - 1]
            temp_output['Evap'] = M.E[1:P.final_day - 1]
            temp_output['Transpiration'] = M.ET_trans[1:P.final_day - 1]
            temp_output['ET'] = M.ET[1:P.final_day - 1]
            temp_output['Pottranspiration'] = M.ET_pot[1:P.final_day - 1]
            temp_output['CumulativePotentialET'] = M.ETcum[1:P.final_day - 1]
            temp_output['CumulativeactualET'] = M.ETcum[1:P.final_day - 1]
            temp_output['ReferenceET'] = M.ET_0[1:P.final_day - 1]
            temp_output['Ky'] = M.Ky[1:P.final_day - 1]
            temp_output['Depletion_total'] = M.Depletion_total[1:P.final_day - 1]
            temp_output['Percent_depletion_total'] = M.Percent_depletion_total[1:P.final_day - 1]
            temp_output['Irrigation'] = M.Irrigation[1:P.final_day - 1]
            temp_output['P'] = M.Ps[1:P.final_day - 1]
            temp_output['Rainfall'] = M.rainfall[1:P.final_day - 1]
            temp_output['Rain_infilt'] = M.Rain_Infilt[1:P.final_day - 1]
            temp_output['Field_capacity_ave'] = M.FCave[1:P.final_day - 1]
            temp_output['Permanent_wilting_point_ave'] = M.PWPave[1:P.final_day - 1]
            temp_output['TAWsum'] = M.TAWsum[1:P.final_day - 1]
            temp_output['RAWsum'] = M.RAWsum[1:P.final_day - 1]
            temp_output['Available_water'] = M.AW[1:P.final_day - 1]
            temp_output['Available_water_depth_ave'] = M.AW[1:P.final_day - 1]
            temp_output['Volumetric_water_content_ave'] = M.VWCave[1:P.final_day - 1]
            temp_output['Total_depth_water'] = M.TDW[1:P.final_day - 1]
            temp_output['Watertableelevation'] = M.TDW[1:P.final_day - 1]
            temp_output['HeightEqu'] = M.Equilibrium_Max[1:P.final_day - 1]
            temp_output['NumberEqu'] = M.zu_Eq_Max[1:P.final_day - 1]
            temp_output['WaterstressKs'] = M.zu_Eq_Max[1:P.final_day - 1]
            temp_output['IrrigationSaltAdded'] = M.Irr_Salt_Added[1:P.final_day - 1]
            temp_output['WasteSaltAdded'] = M.Waste_Salt_Added[1:P.final_day - 1]
            temp_output['SeepageSaltLost'] = M.Seepage[1:P.final_day - 1]
            temp_output['SumOfFluxes'] = M.SumOfFluxes[1:P.final_day - 1]
            temp_output['TotalSalt'] = M.Total_salt[1:P.final_day - 1]
            temp_output['DifferenceInSalt'] = M.SaltDiff[1:P.final_day - 1]
            temp_output['SaltStress'] = M.Ks_salt[1:P.final_day - 1]
            temp_output['NitStress'] = M.Ks_nit[1:P.final_day - 1]
            temp_output['TotalIrrN'] = M.TotalIrrN[1:P.final_day - 1]
            temp_output['TotalDrainN'] = M.TotalDrainN[1:P.final_day - 1]
            temp_output['SumReactionsN'] = M.SumReactionsN[1:P.final_day - 1]
            temp_output['SumFluxReactionsN'] = M.SumFluxReactionsN[1:P.final_day - 1]
            temp_output['TotalMin'] = M.Total_Min[1:P.final_day - 1]
            temp_output['TotalDen'] = M.Total_Den[1:P.final_day - 1]
            temp_output['TotalFer'] = M.Total_Fer[1:P.final_day - 1]
            temp_output['TotalUpt'] = M.Total_Upt[1:P.final_day - 1]
            temp_output['UptakeReq'] = M.UptakeReq[1:P.final_day - 1]
            temp_output['TotalNit'] = M.Total_nit[1:P.final_day - 1]
            temp_output['TotalNitAccum'] = M.Total_nit_accum[1:P.final_day - 1]
            temp_output['NitAve'] = M.N_ave[1:P.final_day - 1]
            temp_output['NOptimalLow'] = M.N_optimal_low[1:P.final_day - 1]
            temp_output['NOptimalHigh'] = M.N_optimal_high[1:P.final_day - 1]
            temp_output['CumMin'] = M.Cum_Min[1:P.final_day - 1]
            temp_output['CumDen'] = M.Cum_Den[1:P.final_day - 1]
            temp_output['CumFer'] = M.Cum_Fer[1:P.final_day - 1]
            temp_output['CumUpt'] = M.Cum_Upt[1:P.final_day - 1]
            temp_output['CumIrrN'] = M.CumIrrN[1:P.final_day - 1]
            temp_output['CumDrnN'] = M.CumDrnN[1:P.final_day - 1]
            temp_output['CumChangeN'] = M.CumChangeN[1:P.final_day - 1]
            
            Number_layers_out = (P.final_day - 1) * (M.Num_layers + 1) + 1

            M.Create_output_layer_arrays(Number_layers_out, P.final_day, M.Num_layers)
            M.date_layer_out[0] = P.DateP[i]
            M.pid_layer_out[0] = P.PlantingIDP[i]
            M.planting_name_layer_out[0] = P.PlantingName[i]
            M.location_layer_out[0] = locate

            temp_layer_output = pd.DataFrame(columns = ['Date', 'Days_after_planting', 'pid', 'plantname', 'Location', 'Layer', 'WaterContent', 'Depletion', 'ActFrac', 'Infiltration',
                                                       'PercentDepletion', 'EC', 'ECe', 'MassSalt', 'Nitrate', 'Nitratemgkg'])
        
            for iii in range(1, P.final_day - 1):
                M.date_layer_out[iii] = dt.fromordinal(M.date_layer_out[iii].toordinal())
            temp_layer_output['Date']=M.date_layer_out[:]
            temp_layer_output['Days_after_planting']=M.Days_after_planting_layer[:]
            temp_layer_output['pid']=M.pid_layer_out[:]
            temp_layer_output['plantname']=M.planting_name_layer_out[:]
            temp_layer_output['Location']=M.location_layer_out[:]
            temp_layer_output['Layer']=M.layer_out[:]
            temp_layer_output['WaterContent']=M.WC_layer_out[:]
            temp_layer_output['Depletion']=M.Depletion_layer_out[:]
            temp_layer_output['ActFrac']=M.ActFrac_layer_out[:]
            temp_layer_output['Infiltration']=M.Infiltration_layer_out[:]
            temp_layer_output['PercentDepletion']=M.Percent_depletion_layer_out[:]
            temp_layer_output['EC']=M.EC_layer_out[:]
            temp_layer_output['ECe']=M.ECe_layer_out[:]
            temp_layer_output['MassSalt']=M.Mass_salt_layer_out[:]
            temp_layer_output['Nitrate']=M.N_layer_out[:]
            temp_layer_output['Nitratemgkg']=M.N_soil_layer_out[:]

            excel_output_path2 = pathprefix + 'output_test3.xlsx'
            writer = pd.ExcelWriter('output_test3.xlsx')
            temp_output.to_excel(writer, 'Output2')
            temp_layer_output.to_excel(writer, 'Output_layers2')
            writer.save()
#            temp_output.to_csv('PythonExport3.csv',sep=',')

            plt.figure(1)
            plt.plot(M.WC)
            plt.ylabel('ml/ml')
            plt.xlabel('DOY')
            plt.show()

            #plt.figure(2)
            #plt.plot(M.
            
            # plt.figure(2)
            # plt.plot(M.rainfall)
            # plt.plot(M.Rain_Infilt)
            # plt.ylabel('mm')
            # plt.xlabel('DOY')
            # plt.show()

            
           
            # temp_layer_output.to_sql('TempLayerOutput', db, if_exists='append', index = 0)
            # temp_output.to_sql('TempOutput',db,if_exists='append', index = 0)

