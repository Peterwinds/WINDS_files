# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 13:38:06 2017

@author: Pete
"""
import pandas as pd
import numpy as np
from datetime import datetime as dt
from datetime import timedelta






def hg(a, b, C, z):
    value = 1
    hg = value
    for i in range(1, 40):
        value = value * (a + i - 1) * (b + i - 1) / (C + i - 1) / i * z
        hg = hg + value
    return hg

def d_total(ts, tr, zt, zwt, alpha, N):
    W = -(np.power(alpha * (zt - zwt) * 100, N))
    zwt = zwt
    hyper = hg(1, 1 - 1 / N, 1 + 1 / N, W / (W - 1))
    d_total = ts * zwt + tr * (zt - zwt) + (ts - tr) * (zt - zwt) * np.power((1 - W), (1 / N - 1)) * hyper
    return d_total

def d_cell_WT(ts, tr, zu, zl, zwt, alpha, N):
    W = -(np.power(alpha * (zu - zwt) * 100, N))
    hyper = hg(1, 1 - 1 / N, 1 + 1 / N, W / (W - 1))
    d_cell_WT = ts * (zwt - zl) + tr * (zu - zwt) + (ts - tr) * (zu - zwt) * np.power((1 - W), (1 / N - 1)) * hyper
    return d_cell_WT

def d_cell_Eq(ts, tr, zu, zl, zwt, alpha, N):
    WU = -(np.power(alpha * (zu - zwt) * 100, N))
    wl = -(np.power(alpha * (zl - zwt) * 100, N))
    hyperu = hg(1, 1 - 1 / N, 1 + 1 / N, WU / (WU - 1))
    hyperl = hg(1, 1 - 1 / N, 1 + 1 / N, wl / (wl - 1))
    d_cell_Eq = tr * (zu - zl) + (ts - tr) * (zu - zwt) * np.power((1 - WU), (1 / N - 1)) * hyperu - (ts - tr) * (zl - zwt) * np.power((1 - wl), (1 / N - 1)) * hyperl
    return d_cell_Eq
    
# The following function, if true, reads the weather from the SQL database, otherwise from Active_year_weather in the spreadsheet
# The data is read into a panda array called WeatherArray and returned
    
 
class weather:
    def __init__(self, WeatherArray):
        self.NumDays = self.EndDOY - self.SimStartDOY 
        self.w_date = WeatherArray['Date']                   #the date of the record, mm/dd/yyyy
        self.RH_min = np.array(WeatherArray['rh_min'])                      #Minimum relative humidity, percent, decimal, 8 digits
        self.WS = np.array(WeatherArray['WS'])                  #wind speed, m/sec, decimal, 8 digits
        self.ET_0 = np.array(WeatherArray['ETo (mm)'])                        #Reference ET, mm/day, decimal, 8 digits
        self.rainfall = np.array(WeatherArray['Rain (mm)'])              #Rainfall, mm/day, decimal, 8 digits
        self.Tmax = np.array(WeatherArray['tmax'])                      #Maximum temperature, Celsius, decimal, 8 digits
        self.Tmin = np.array(WeatherArray['tmin'])                      #Minimum temperature, Celsius, decimal, 8 digits
        self.Rain_time = np.array(WeatherArray['Rain time (hr)'])                 #Rain time (time of storm), hr, decimal, 5 digits

class RS_daily:
    def __init__(self, RS_dailyArray):
        self.NDVI_DOY = np.array(RS_dailyArray['NDVI_DOY'])           #The name of the planting (text), normally includes owner name, field name, crop, and year, such as GaryField1Wheat2017
        self.NDVI_daily = np.array(RS_dailyArray['NDVI'])           #The name of the planting (text), normally includes owner name, field name, crop, and year, such as GaryField1Wheat2017

class ET_daily:
    def __init__(self, ET_dailyArray):
        self.Daily_DOY = np.array(ET_dailyArray['DOY'])                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
        self.Daily_Kcb = np.array(ET_dailyArray['Kcb'])
        self.Daily_one_minus_c = np.array(ET_dailyArray['one_minus_c'])
        self.Daily_Crop_height = np.array(ET_dailyArray['Crop_height'])
        self.Daily_Root = np.array(ET_dailyArray['Root'])

class plantings:
    def __init__(self, PlantingArray):
        self.PlantingNum = PlantingArray['Planting num']                             #The long random identifier that represents the planting, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.PlantingName = PlantingArray['Planting name']                #The name of the planting (text), normally includes owner name, field name, crop, and year, such as GaryField1Wheat2017
        self.RunPlanting = PlantingArray['Run planting']                 #A Boolean that specifies whether a simulation should take place
        self.PlantingYear = PlantingArray['Planting year']                 #A Boolean that specifies whether a simulation should take place
        self.PlantingDOY = PlantingArray['Planting DOY']
        self.PlantingDate = PlantingArray['Planting date']
        self.StartDOY = PlantingArray['DOY start']                             #The long random identifier that represents the planting, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.StartDate = PlantingArray['Date start']                #The name of the planting (text), normally includes owner name, field name, crop, and year, such as GaryField1Wheat2017
        self.EndDOY = PlantingArray['DOY end']
        self.EndDate = PlantingArray['Date end']                   #The number of days ahead of the present day (integer) for which the model should predict growth, water use, and irrigations
        self.FieldNum = PlantingArray['Field num']                                #The long random identifier that represents the field, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.FieldName = PlantingArray['Field name']                             #The name of the field (text), which normally includes the owner and field, such as GaryField1
        self.AccountNum = PlantingArray['Account name']                              #The long random identifier that represents the owner, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.AccountName = PlantingArray['Account name']        #The number of locations (integer) that are modeled in an individual planting
        self.CropNum = PlantingArray['Crop num']                         #The name of the crop (text string), which initializes the following parameters, for example, GuayuleAZ2Spring
        self.CropName = PlantingArray['Crop name']                              #The crop ID, which initializes the following parameters, but the parameters can be subsequently changed, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9 
        self.WeatherNum = PlantingArray['Weather num']       #The name of the weather station (text) to be used for this simulation, for example, Maricopa represents the AZMET Maricopa station
        self.WeatherName = PlantingArray['Weather name']          #the station ID of the record, mm/dd/yyyy
        self.WeatherSheetName = PlantingArray['Weather sheet name']
        self.ETFractionsNum = PlantingArray['ET fractions num']                              #The crop ID, which initializes the following parameters, but the parameters can be subsequently changed, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9 
        self.WettingFractionsNum = PlantingArray['Wetting fractions num']
        self.SimStartDOY = PlantingArray['SimStartDOY']
        self.SimStartDate = PlantingArray['SimStartDate']
        self.IrrNum = PlantingArray['Irr num']
        self.IrrName = PlantingArray['Irr name']
        self.IrrSection = PlantingArray['Irr section']
        self.IrrFirstRow = PlantingArray['Irr first row']
        self.NumIrrigations = PlantingArray['Num irrigations']
        self.NumNeutron = PlantingArray['Num neutron']
        self.FirstNeutronRow = PlantingArray['First neutron row']
        self.IP = PlantingArray['Initial_days']                                   #The number of days after planting (integer) that the development phase (IncreasingET) begins
        self.DP = PlantingArray['Development_days']                     #The number of days after planting (integer) that the Midseason (maxET) phase begins
        self.MP = PlantingArray['Midseason_days']                 #The number of days after planting (integer) that the Late Season decline phase begins
        self.LP = PlantingArray['Lateseason_days']               #The number of days after planting (integer) that the End Season (low ET) phase begins. Some crops so not have this phase
        self.EP = PlantingArray['Ending_season_days']               #The number of days after planting (integer) that the End Season (low ET) phase begins. Some crops so not have this phase
        self.DP_2 = PlantingArray['Development_days_year_2']               #The number of days after planting (integer) that the End Season (low ET) phase begins. Some crops so not have this phase
        self.MP_2 = PlantingArray['Midseason_days_year_2']               #The number of days after planting (integer) that the End Season (low ET) phase begins. Some crops so not have this phase
        self.LP_2 = PlantingArray['Lateseason_days_year_2']               #The number of days after planting (integer) that the End Season (low ET) phase begins. Some crops so not have this phase
        self.CropEndDOY = PlantingArray['Crop_end_DOY']                         #The number of days after planting (integer) that the harvest takes place.
        self.IKcb = PlantingArray['Initial_Kcb']                     #The crop coefficient during the initial stage (decimal, 5 digits), based on dual crop coefficient and FAO56
        self.MKcb = PlantingArray['Midseason_Kcb']                             #The crop coefficient during the midseason stage (decimal, 5 digits), based on dual crop coefficient and FAO56. The model increases the crop coefficient during the development stage
        self.EKcb = PlantingArray['Endseason_ Kcb']                             #The crop coefficient during the midseason stage (decimal, 5 digits), based on dual crop coefficient and FAO56. The model increases the crop coefficient during the development stage
        self.MKcb_2 = PlantingArray['Midseason_Kcb_year_2']                             #The crop coefficient during the endseason stage (decimal, 5 digits), based on dual crop coefficient and FAO56. The model decreases the crop coefficient during the lateseason stage
        self.EKcb_2 = PlantingArray['Endseason_Kcb_year_2']                             #The crop coefficient during the endseason stage (decimal, 5 digits), based on dual crop coefficient and FAO56. The model decreases the crop coefficient during the lateseason stage
        self.ICP = PlantingArray['Initial_1-C_root_height_Days']                                   #The number of days after planting (integer) that the development phase (IncreasingET) begins
        self.DCP = PlantingArray['Development_1-C_root_height_Days']                     #The number of days after planting (integer) that the Midseason (maxET) phase begins
        self.MCP = PlantingArray['Midseason_1-C_root_height_Days']                     #The number of days after planting (integer) that the Midseason (maxET) phase begins
        self.ECP = PlantingArray['Endseason_1-C_Days']                 #The number of days after planting (integer) that the Late Season decline phase begins
        self.IC = PlantingArray['Initial_1-C']                     #The number of days after planting (integer) that the Midseason (maxET) phase begins
        self.MC = PlantingArray['Midseason_1-C']                 #The number of days after planting (integer) that the Late Season decline phase begins
        self.EC = PlantingArray['Endseason_1-C']                 #The number of days after planting (integer) that the Late Season decline phase begins
        self.ERH = PlantingArray['Endseason_root_and_height_days']                 #The number of days after planting (integer) that the Late Season decline phase begins
        self.IR = PlantingArray['Initial_root_z']                     #The number of days after planting (integer) that the Midseason (maxET) phase begins
        self.DR = PlantingArray['Development_root_z']                 #The number of days after planting (integer) that the Late Season decline phase begins
        self.FR = PlantingArray['Final_root_z']                 #The number of days after planting (integer) that the Late Season decline phase begins
        self.minH = PlantingArray['Minimum_plant_height']                   #The initial crop height in the model (decimal, 5 digits, meters).
        self.maxH = PlantingArray['Maximum_plant_height']                   #The final crop height (decimal, 5 digits, meters). The model grows the crop during the development phase
        self.midH = PlantingArray['Midseason_plant_height']                   #The initial crop height in the model (decimal, 5 digits, meters).
        self.endH = PlantingArray['Endseason_plant_height']                   #The final crop height (decimal, 5 digits, meters). The model grows the crop during the development phase
        self.calculate_one_minus_c = PlantingArray['calculate_one_minus_c']                   #The final crop height (decimal, 5 digits, meters). The model grows the crop during the development phase
        self.Select_Kcb = PlantingArray['Kc_calculation_procedure']                          #The management allowable depletion (MAD) value from table 22 in FAO56. It is then adjusted for evaporation rate. fraction, 5 digits
        self.Irrigation_partition_1 = PlantingArray['Irrigation_partition_1']
        self.Irrigation_partition_2 = PlantingArray['Irrigation_partition_2']
        self.Irrigation_partition_3 = PlantingArray['Irrigation_partition_3']
        self.Irrigation_partition_4 = PlantingArray['Irrigation_partition_4']
        self.Nitrogen_simulation = PlantingArray['Simulate_nitrogen']
        self.Salinity_simulation = PlantingArray['Simulate_salinity']
        self.Rainfall_simulation = PlantingArray['Simulate_rainfall']
        self.type_irrigation = PlantingArray['Type_irrigation']
        self.Rain_infiltration_calculation = PlantingArray['Rain_infiltration_calculation']    #Boolean to determine if infiltration is calculated with Green ampt or SCS or is just assigned the rain_partition as the fraction

        self.Neglect_upper_layer_in_depletion_calc = PlantingArray['NeglectEvaporationLayerDepletion'] #Boolean that is true if upper layer is  not included in the calculation of percent depletion
        self.pTable22 = PlantingArray['MAD_p_Table22']                          #The management allowable depletion (MAD) value from table 22 in FAO56. It is then adjusted for evaporation rate. fraction, 5 digits
        self.adjust_P = PlantingArray['Adjust_p']
        self.Max_yield = PlantingArray['Max_yield']                          #in case of economic forecast, this is in the model, but it is not currently used. (kg/ha) (Real number, 8 digits)
        self.Ky_initial = PlantingArray['Initial_crop_sensitivity']             #The sensitivity of the crop yield to water stress during initial phase (Units can vary dramatically for different scenarios) (Decimal number, 5 digits)
        self.Ky_DP = PlantingArray['Development_crop_sensitivity']                  #The sensitivity of the crop yield to water stress during development phase (Units can vary dramatically for different scenarios) (Decimal number, 5 digits)
        self.Ky_MP = PlantingArray['Midseason_crop_sensitivity']                  #The sensitivity of the crop yield to water stress during midseason phase (Units can vary dramatically for different scenarios) (Decimal number, 5 digits)
        self.Ky_late = PlantingArray['Lateseason_crop_sensitivity']                #The sensitivity of the crop yield to water stress during lateseason phase (Units can vary dramatically for different scenarios) (Decimal number, 5 digits)
        self.Ky_end = PlantingArray['Endseason_crop_sensitivity']                  #The sensitivity of the crop yield to water stress during endseason phase (Units can vary dramatically for different scenarios) (Decimal number, 5 digits)
        self.Seal_bottom = PlantingArray['seal_bottom']                     #Boolean that is true for some class assignments in which soil physics calculations are evaluated
        self.No_infiltration = PlantingArray['no_infiltration']             #Boolean that is true for some class assignments in which soil physics calculations are evaluated
        self.No_ET = PlantingArray['no_et']                                 #Boolean that is true for some class assignments in which soil physics calculations are evaluated
        self.No_stress_reduction =PlantingArray['no_stress_reduction']                #Boolean that is true if water or salt stress is not allowed to reduce plant uptake of water .
        self.No_ET_frac_Adjustment = PlantingArray['NoETFracAdjustment']   #Boolean that is true if there is no adjustment of the layers from which water is extracted due to low water content in some layers. 
        self.No_Redistribution = PlantingArray['NoRedistributionRichardsEqn']  #Boolean that is true if Richard's equation dateddddds not redistribute water between irrigations
        self.Field_capacity_restriction = PlantingArray['Field_capacity_restriction']   #Boolean that is true if water in a layer is allowed to be greater than field capacity during infiltration
        self.Irrigation_efficiency = PlantingArray['Irrigation_efficiency']
        self.LeachingFraction = PlantingArray['Leaching_fraction']          #The fraction of water that leaches below the root zone when calculating the next irrigation depth. (decimal, 5 digits)      
        self.Eliminate_surface_evaporation = PlantingArray['Eliminate_surface_evaporation']
        self.DOY_to_eliminate_surface_evaporation = PlantingArray['DOY_to_eliminate_surface_evaporation']
        self.Num_wetting_phases = PlantingArray['Num_wetting_phases']  #this is the number of wetting phases, such as 
        self.FW_phase_DOY_1 = PlantingArray['Wetting_switch_day_1']
        self.FW_phase_DOY_2 = PlantingArray['Wetting_switch_day_2']
        self.FW_phase_DOY_3 = PlantingArray['Wetting_switch_day_3']
        self.FW_phase_DOY_4 = PlantingArray['Wetting_switch_day_4']
        self.FW_phase_DOY_5 = PlantingArray['Wetting_switch_day_5']
        self.FW_Typ_1 = PlantingArray['Irrigation_type_1']
        self.FW_Type_2 = PlantingArray['Irrigation_type_2']
        self.FW_Type_3 = PlantingArray['Irrigation_type_3']
        self.FW_Type_4 = PlantingArray['Irrigation_type_4']
        self.FW_Type_5 = PlantingArray['Irrigation_type_5']

#Salinity parameters
        self.ECiw = PlantingArray['Irrigation_water_EC'] #Salinity of the irrigation water, ranges from 0 to 5, decimal, dS/m, 5 digits
        self.Max_S = PlantingArray['Maximum_soluble_EC_act'] #Maximum solubility of salts in soil solution, actual salinity, EC, dS/m, 8 digits
        self.Waste_app_DOY = PlantingArray['DOY of salinity application'] #Date of waste (manure) application   mm/dd/yyyy
        self.Waste_sal_kgha = PlantingArray['Rate_of_salinity_application_kgha'] #The kg/ha of salts in waste applied to the field, decimal, 5 digits
        self.Waste_dissolution = PlantingArray['Rate_of_dissolution'] #1/number of days for waste to dissolve, fraction, 5 digits
        self.ECe_thresh = PlantingArray['ECe_thresh'] #the threshhold ECe, below which there is no reduction in growth, real number, 5 digits
        self.b_sal = PlantingArray['b_sal'] #The slope of the yield, ECe line
        self.WasteAppTF = PlantingArray['WasteAppTF']
        
#Nitrogen parameters
        self.nit_dissolution = PlantingArray['Fertilizer_dissolution_rate'] #1/days to dissolve into soil, fraction, 5 digits
        self.Kmnl = PlantingArray['Kmnl'] # Mineralization constant /day, fraction 8 digits
        self.qlow = PlantingArray['qlow'] # Lower threshold water content for mineralization, ml/ml, fraction 5 digits
        self.qhigh = PlantingArray['qhigh'] # Upper threshold water content for mineralization, ml/ml, fraction 5 digits
        self.Qtemp = PlantingArray['Qtemp'] # Rate of change associated with 10 degree C change in temperature, dimensionless, decimal, 8 digits
        self.Onmax = PlantingArray['Surface_organic_matter'] # Organic matter content in top layer of soil, micrograms/gram of soil, real number, 8 digits
        self.alpha = PlantingArray['alpha'] # Constant for decrease in organic matter content with depth, dimensionless, fraction 8 digits
        self.Kden = PlantingArray['Kden'] #Denitrification constant, /day, fraction, 8 digits
        self.Nmin = PlantingArray['Nmin'] # M-Menton coefficient, (mg NO3-N)/(kg soil), decimal, 5 digits
        self.Km = PlantingArray['Km'] #M-Menton coefficient, (mg NO3-N)/(kg soil), decimal, 5 digits
        self.Frac_NO3 = PlantingArray['Fraction_N_req_taken_as_nitrate'] # Fraction of total nitrogen req. taken up as nitrate (vs. ammonium), decimal, 5 digits
        self.Seasonal_N_uptake = PlantingArray['Seas_N_requirement'] # Total seasonal nitrogen requirement, kg/ha, real, 5 digits
        self.Fert1_DOY = PlantingArray['DOY_of_Fert_app_1'] #Date on which first fertilization takes place mo/da/year
        self.Fert2_DOY = PlantingArray['DOY_of_Fert_app_2'] #Date on which second fertilization takes place  mo/da/year
        self.Fert3_DOY = PlantingArray['DOY_of_Fert_app_3'] #Date on which third fertilization takes place mo/da/year
        self.Fert1_rate = PlantingArray['Rate_of_application_1'] #Rate of fertilizer application 1 (kg/ha) real, 8 digits
        self.Fert2_rate = PlantingArray['Rate_of_application_2'] #Rate of fertilizer application 2 (kg/ha) real, 8 digits
        self.Fert3_rate = PlantingArray['Rate_of_application_3'] #Rate of fertilizer application 3 (kg/ha) real, 8 digits
        self.Fert_depth = PlantingArray['Depth_of_fertilizer_application']
        self.Number_fertilizations = PlantingArray['Number_of_fertilization_events'] #Number of times fertilizer is applied, Integer
        self.Constant_Km = 1
        self.Constant_Km_1 = PlantingArray['Read_in_daily_Km_values']
        self.Constant_Km_2 = PlantingArray['Adjust_Km_in_program_for_opt']
        self.Constant_Km_3 = PlantingArray['Use_constant_Km_values']
        self.N_soil_optimal = PlantingArray['Optimal_soil_nitrate_concentration'] #Optimal soil nitrate concentration for maximum yield, real, 5 digits
        self.Frac_greater_saturation = PlantingArray['Sat_uptake_above_optimal'] #Fraction of nitrate uptake greater than optimal at saturated nitrate condition, fraction, 5 digits
        self.Upper_temp_adjust = PlantingArray['Upper_limit_of_temperature_adjust'] #Upper limit of of min den adjustment for temperature, fraction, 5 digits
        self.Upper_temp_reduction = PlantingArray['Rate_of_decrease_after_upper_lim'] #Reduction rate of temp adjustment above upper limit, fraction, 5 digits
        self.Niw = PlantingArray['Irrigation_nitrate_conc'] #Irrigation water nitrate concentration (constant for season) mg/L
        self.Range_frac = PlantingArray['Fraction_of_optimal_range'] #This is used after M-Menton to assign yield decrease based on change from optimal, fraction, 5 digits
        self.Yield_decrease_below = PlantingArray['Rate_of_yield_decrease_below_N'] #Rate of yield decrease below optimal N uptake range
        self.Yield_decrease_above = PlantingArray['Rate_of_yield_decrease_above_N'] #Rate of yield decrease above optimal N uptake range
        self.Dont_fertilize_evap_layer = PlantingArray['Dont_fertilize_evap_layer']
        self.tm = PlantingArray['tm'] #Rate of yield decrease above optimal N uptake range
        self.T_bar = PlantingArray['T_bar'] #Rate of yield decrease above optimal N uptake range
        self.A0 = PlantingArray['A0'] #Rate of yield decrease above optimal N uptake range
        self.diff = PlantingArray['Soil_thermal_diffusivity']

class fields:
    def __init__(self, Field_Array):
        print(Field_Array)
        self.FieldNum = Field_Array['Field num']                               #The long random identifier that represents the field, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.FieldName = Field_Array['Field name']                             #The name of the field (text), which normally includes the owner and field, such as GaryField1
        self.AccountNum = int(Field_Array['Account num'])                       #The long random identifier that represents the owner, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.AccountName = Field_Array['Account name']                              #The long random identifier that represents the owner, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.Num_layers = int(Field_Array['Num_layers'])                        #The number of layers (integer), not including the evaporation layer, that will be modeled in the current planting, which may be less than the number of field layers
        self.SAV = float(Field_Array['SAV'])                                    #Suction at wetting front (cm) decimal 5 digits, used for Green Ampt infiltration
        self.H0 = float(Field_Array['H0'])                   #Depth of water ponded in field if Green Ampt is used to calculate infiltration (I do not know why this is cm but it is) (decimal 5 digits)
        self.REW = float(Field_Array['REW']) #Readily evaporable water in the evaporation layer (read in as cm, decimal, 8 digits)
        self.TEW = float(Field_Array['TEW']) #Total evaporable water in the evaporation layer (read in as cm, decimal 8 digits)
        self.NRCS_a = float(Field_Array['NRCS_a'])                                  #The SCS a value for Kostiakov infiltration (decimal, 5 digits)
        self.NRCS_b = float(Field_Array['NRCS_b'])                                  #The SCS b value for Kostiakov infiltration (decimal, 5 digits)
        self.Intake_family = float(Field_Array['Intake_family'])              #Boolean to determine whether SCS curve numbers are used to determine infiltration
        self.Water_table_simulation = int(Field_Array['Water_table_simulation']) #Boolean that is true if there is a water table and leaching is restricted by drainage system
        self.Begin_equil = int(Field_Array['Begin_with_soil_water_in_equilibrium_with_WT'])
        self.Init_zwt = float(Field_Array['Initial_water_table_elevation']) #Initial height of the water table on the day before planting (decimal, 5 digits, meters)
        self.Drain_elevation = float(Field_Array['Drain_Elevation']) #Elevation of the drain within  the soil profile (decimal, 5 digits, meters)
        self.Kirkham_rate = int(Field_Array['Drainage_controlled_by_Kirkham']) #Boolean that is true if the Kirkham algorithm is used to calculate drainage rate
        self.Drain_rate = float(Field_Array['Multiplier_for_linear']) #A multiplier that calculates drainage rate as a function of water table height (decimal, 8 digits, dimensionless) if not using Kirkham method
        self.Drain_imp = float(Field_Array['Drain_elevation_above_impermeable']) #The elevation of the drain above the impermeable layer (decimal, 8 digits, meters), currently used in Kirkham form
        self.L_drain = float(Field_Array['Distance_between_drains']) #The distance between drains (decimal, 8 digits, m), currently used in Kirkham form
        self.Keff_horizontal = float(Field_Array['Effective_lateral_K']) #Effective horizontal conductivity (m/day, 8 digits, decimal), currently used in Kirkham form
        self.re = float(Field_Array['Effective_drain_radius']) #Effective drain radius, m (decimal, 8 digits)
        self.Horizontal_distance_from_drain = float(Field_Array['Horizontal_distance_from_drain'])
        self.Kirkham_F = float(Field_Array['Kirkham_f']) #Calculated value based on drainage geometry (m/day, decimal, 8 digits)
        self.Max_layer_of_equilibrium_layers = float(Field_Array['Max_equilibrium_layer'])
        self.Equil_max_init = float(Field_Array['Starting_equilibrium_layer']) #this is the fraction of water in the soil, above which a layer is considered in equilibrium with the water table, fraction, 5 digits
        self.Keep_equilibrium_layers_out_of_root_zone = float(Field_Array['Keep_equilibrium_layers_out_of_root_zone']) #Boolean that is true if equilibrium layers in drainage simulation are not allowed to enter the root zone
        self.Fraction_of_saturation_for_equilibrium = float(Field_Array['Fraction_of_saturation_for_equilibrium']) #Fraction of saturation for equilibrium, fraction, 5 digits
        self.Continue_drainage_rate = float(Field_Array['Continue_drainage_rate'])
        self.Rain_partition = 1
        
    
class soil:
    def __init__(self, SoilArray):
        SoilArray = SoilArray.sort_values('Layer')
        self.FieldNum = np.array(SoilArray['Field num'])                            #The long random identifier that represents the field, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.LayerNumber = np.array(SoilArray['Layer'])
        self.Depth = np.array(SoilArray['Depth'])  
        self.InitWC = SoilArray['InitWC']
        self.switch_wetting_phase = 400
        self.FC = np.array(SoilArray['FC'])                                #The field capactiy of the layer, percentage, decimal, 8 digits
        self.PWP = np.array(SoilArray['PWP'])                              #The permanent wilting point of the layer, percentage, decimal, 8 digits
        self.ts = np.array(SoilArray['Sat'])                               #The saturated water content of the layer, percentage, decimal, 8 digits
        self.Ksat = np.array(SoilArray['Ksat'])                            #The saturated vertical hydraulic conductivity of the layer, percentage, decimal, 8 digits
        self.av = np.array(SoilArray['alphav'])                            #Van Genuchten alpha, decimal, 8 digits
        self.nv = np.array(SoilArray['nv'])                                #Van Genuchten n, decimal, 8 digits
        self.tr = np.array(SoilArray['ResWC'])                             #Residual water content, from Brooks Corey, decimal, 8 digits
        self.Lv = np.array(SoilArray['Lv'])                                #Van Genuchten L, decimal, 8 digits
        self.Ko = np.array(SoilArray['Ko'])                                #Matching Ko (conductivity) used in Van Genuchten hydraulic conductivity equation, decimal, 8 digits
        self.InitECe = np.array(SoilArray['ECe_init'])                     #The initial ECe at the beginning of the season in the layer, dS/m, decimal, 8 digits
        self.InitsoilN = np.array(SoilArray['Initial_soil_N'])              #Initial nitrogen concentration mg/L soil
        self.Org = np.array(SoilArray['Organic_matter'])              #Initial nitrogen concentration mg/L soil
        self.dz = np.zeros([len(self.Depth)+1])                               #the thickness of the layer, m, decimal, 8 digits
        self.Z_top_ave = np.zeros([len(self.Depth)+1])                        #Half of the thickness of two layers, m, decimal, 8 digits
        self.Total_TAW = np.zeros([len(self.Depth)+1])                        #Total available water from the top of profile down to and including the layer, m, decimal, 8 digits
        self.AWC = np.zeros([len(self.Depth)+1])                              #Readily available water from the top of profile down to and including the layer, m, decimal, 8 digits
        self.RAW = np.zeros([len(self.Depth)+1])                              #Readily available water from the top of profile down to and including the layer, m, decimal, 8 digits
        self.TAW = np.zeros([len(self.Depth)+1])                              #Readily available water from the top of profile down to and including the layer, m, decimal, 8 digits
        self.bd = np.zeros([len(self.Depth)+1])                               #Bulk density of the layer, g/cm3, decimal, 8 digits
        self.zu = np.zeros([len(self.Depth)+1])
        self.zl = np.zeros([len(self.Depth)+1])
        self.Zave = np.zeros([len(self.Depth)+1])
        self.Fw = np.zeros([len(self.Depth)+1])
        self.Fw_y = np.zeros([len(self.Depth)+1])
        self.Active_EC = np.zeros([len(self.Depth)+1])
        self.Active_n = np.zeros([len(self.Depth)+1])
        self.Active_WC = np.zeros([len(self.Depth)+1])
        self.Active_dry_n = np.zeros([len(self.Depth)+1])
        self.mv = np.zeros([len(self.Depth)+1])
        
        self.InitWC = self.InitWC / 100
        self.FC = self.FC / 100
        self.PWP = self.PWP/ 100
        self.AWC = self.FC - self.PWP
        self.ts = self.ts / 100
        self.bd = (1-self.ts)*2.65
        self.tr = self.tr / 100
        self.Ksat = self.Ksat / 100
        self.Ko = self.Ko / 100
        self.Depth = self.Depth / 100
        self.dz[0] = self.Depth[0]
        self.Z_top_ave[0] = self.Depth[0] / 2
        print(self.Num_layers, self.Depth, 'here')
        for i in range(1, self.Num_layers+1):
            self.dz[i] = self.Depth[i] - self.Depth[i+1]
        self.dz[self.Num_layers + 1] = self.Depth[self.Num_layers + 1]
        for i in range(1, self.Num_layers + 2):
            self.Z_top_ave[i] = self.Depth[i] - self.dz[i]/2
        for i in range(1, self.Num_layers + 2):
            self.TAW[i] = self.AWC[i] * self.dz[i]
            self.RAW[i] = self.TAW[i] * self.pTable22
            self.mv[i] = 1 - 1 / self.nv[i]
        self.Total_TAW[0] = 0
        for i in range(1, self.Num_layers + 2):
            self.Total_TAW[i] = self.Total_TAW[i - 1] + self.TAW[i]
        for i in range(1, self.Num_layers + 2):
               if (i == 1):
                    self.zu[i] = self.dz[1]
                    self.zl[i] = 0
               elif i == self.Num_layers + 1:
                    self.zu[i] = self.Depth[1]
                    self.zl[i] = self.Depth[1] - self.Depth[i]
               else:
                    self.zu[i] = self.Depth[1] - self.Depth[i+1]
                    self.zl[i] = self.Depth[1] - self.Depth[i]
               self.Zave[i] = (self.zu[i] + self.zl[i]) / 2 
       

class ET_fractions:
    def __init__(self, ET_Fractions_Array):

        self.ET_fractions = np.zeros((self.Num_layers + 2, self.Num_layers + 2))

        for j in range(1, int(self.Num_layers) + 2):
            for k in range(1, int(self.Num_layers) + 2):
                    self.ET_fractions[j][k] = ET_Fractions_Array.iloc[13-k][self.Num_layers + 5-j]
                    
                    
class Wetting_fractions:
    def __init__(self, Wetting_Fractions_Array):

        self.FW_layers = np.zeros((self.EndDOY + 1, self.Num_layers + 2))

        for j in range(1, self.EndDOY + 1):
            for k in range(1, int(self.Num_layers) + 2):
                if j < self.FW_phase_DOY_1:
                    self.FW_layers[j][k] = Wetting_Fractions_Array.iloc[0][4+k]
                elif j < self.FW_phase_DOY2:
                    self.FW_layers[j][k] = Wetting_Fractions_Array.iloc[1][4+k]
                elif j < self.FW_phase_DOY3:
                    self.FW_layers[j][k] = Wetting_Fractions_Array.iloc[2][4+k]
                elif j < self.FW_phase_DOY4:
                    self.FW_layers[j][k] = Wetting_Fractions_Array.iloc[3][4+k]
                else:
                    self.FW_layers[j][k] = Wetting_Fractions_Array.iloc[4][4+k]
                    
class irrigation:
    def __init__(self, IrrigationArray):
        self.Irrigation_depth = np.zeros(self.NumDays+1)
        self.Depth_ref = np.array(IrrigationArray['Ref_mm']) 
        Sec_name = 'Sec_' + str(self.IrrSection)
        self.Multiplier = np.array(IrrigationArray[Sec_name])
        self.DOY = np.array(IrrigationArray['DOY'])
        for i in range(0, len(IrrigationArray)):
            self.Irrigation_depth[self.DOY[i]-self.StartDOY] = self.Depth_ref[i]*self.Multiplier[i]


class status:
    def __init__(self, Status_Array):

        self.Water_Start = np.zeros(len(Status_Array)+1)                              #Readily available water from the top of profile down to and including the layer, m, decimal, 8 digits
        self.Salinity_Start = np.zeros(len(Status_Array)+1)                              #Readily available water from the top of profile down to and including the layer, m, decimal, 8 digits
        self.Nitrogen_Start = np.zeros(len(Status_Array)+1)                              #Readily available water from the top of profile down to and including the layer, m, decimal, 8 digits

        for k in range(1, int(self.Num_layers) + 2):
            self.Water_Start[k] = Status_Array.iloc[k][4]
            self.Salinity_Start[k] = Status_Array.iloc[k][5]
            self.Nitrogen_Start[k] = Status_Array.iloc[k][6]
                    
class output:
    def __init__(self, Output_Array):
        self.datelist = [self.StartDate + timedelta(days = x) for x in range(0, 1 + self.NumDays)]
        self.DOY = np.zeros(self.NumDays)
        self.DOY[0] = self.StartDOY
        self.PlantingIDO = np.empty(self.NumDays, dtype = object)                  #The long random identifier that represents the planting, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.PlantingNameO = np.empty(self.NumDays, dtype = object)                 #The name of the planting (text), normally includes owner name, field name, crop, and year, such as GaryField1Wheat2017
        self.LocationO = np.zeros(self.NumDays)                     #The position (location) in the planting, integer
        self.NDVI = np.zeros(self.NumDays)                          #The NDVI measured by remote sensing for the location, decimal 8 digits
        self.NDVI_regression = np.zeros(self.NumDays)               #The interpolated NDVI by regression, decimal 8 digits
        self.Kcb_NDVI = np.zeros(self.NumDays)                      #The Kcb calculated from the NDVI regression, decimal 8 digits
        self.Kcb_calculated = np.zeros(self.NumDays)                #The Kcb calculated from FAO equations, decimal 8 digits
        self.Kcb = np.zeros(self.NumDays)                           #The Kcb used in the program, decimal 8 digits                            
        self.one_minus_c = np.zeros(self.NumDays)                   #The one minus the fraction of the canopy area (ground area fraction), decimal 8 digits
        self.Crop_height = np.zeros(self.NumDays)        #The h8 of the crop in meters, decimal 8 digits
        self.Root = np.zeros(self.NumDays)                #The depth of the root in meters, decimal 8 digits
        self.Ke = np.zeros(self.NumDays)                       #The evaporation coefficient, decimal 8 digits
        self.E = np.zeros(self.NumDays)                       #The evaporation from the soil surface, mm, decimal 8 digits
        self.ET_trans = np.zeros(self.NumDays)        #The transpiration from the plant, mm, decimal 8 digits
        self.ET = np.zeros(self.NumDays)                         #The evapotranspiration from the crop, mm, decimal 8 digits
        self.ET_pot = np.zeros(self.NumDays)         #The potential transpiration from the crop, decimal 8 digits
        self.ETcum = np.zeros(self.NumDays)      #The cumulative potential evapotranspiration for the season, decimal 8 digits
        self.ET_actual_cum = np.zeros(self.NumDays)#The cumulative actual evapotranspiration for the season, decimal 8 digits
        self.ET_refO = np.zeros(self.NumDays)
        self.Ky = np.zeros(self.NumDays)                           #The daily crop sensitivity to water stress
        self.Depletion_total = np.zeros(self.NumDays)  #The depth of depletion, mm, decimal 8 digits
        self.Percent_depletion_total = np.zeros(self.NumDays) #The percent depletion, percentage, decimal 8 digits
        self.Irrigation = np.zeros(self.NumDays)           #The depth of irrigation applied, mm, decimal 8 digits
        self.P = np.zeros(self.NumDays)                           #The management allowed depletion, fraction, decimal 8 digits
        self.Rainfall = np.zeros(self.NumDays)                #The rainfall depth, mm, decimal 8 digits
        self.Rain_Infilt = np.zeros(self.NumDays)
        self.wet_Rain_Infilt = np.zeros(self.NumDays)
        self.FCave = np.zeros(self.NumDays)                 #The average field capacity for soil profilel, decimal 8 digits
        self.PWPave = np.zeros(self.NumDays)                #The average permanent wilting point for soil profilel, decimal 8 digits
        self.TAWsum = np.zeros(self.NumDays)                #The average permanent wilting point for soil profilel, decimal 8 digits
        self.RAWsum = np.zeros(self.NumDays)                #The average permanent wilting point for soil profilel, decimal 8 digits
        self.AW =np.zeros(self.NumDays)                 #The amount of water left in the soil for the plant, mm, decimal 8 digits
        self.VWCave = np.zeros(self.NumDays)                #The average volumetric water content for the soil profile, fraction, decimal 8 digits
        self.TDW = np.zeros(self.NumDays)                   #The total depth of water in the soil profile, mm, decimal 8 digits           
        self.zwt = np.zeros(self.NumDays)                   #The water table elevation above the datum, m, decimal 8 digits
        self.Equilibrium_Max = np.zeros(self.NumDays)       #The number of the upper layer that is in equilibrium with the water table, integer
        self.zu_Eq_Max = np.zeros(self.NumDays)       #The number of the upper layer that is in equilibrium with the water table, integer
        self.Ks = np.zeros(self.NumDays)                    #Ks which is the fraction of decreased evapotranspiration due to limited water, decimal, 8 digits
        self.Irr_Sal = np.zeros(self.NumDays)              #The salt mass added by irrigation water, mg/L m, decimal, 8 digits
        self.Waste_sal = np.zeros(self.NumDays)            #The salt mass added by waste, mg/L m, decimal, 8 digits
        self.Irr_Salt_Added = np.zeros(self.NumDays)               #The salt mass lost by seepage, mg/L m, decimal, 8 digits
        self.Waste_Salt_Added = np.zeros(self.NumDays)           #The mass balance of all salts added or lost, mg/L m, decimal, 8 digits
        self.Seepage = np.zeros(self.NumDays)               #The salt mass lost by seepage, mg/L m, decimal, 8 digits
        self.SumOfFluxes = np.zeros(self.NumDays)           #The mass balance of all salts added or lost, mg/L m, decimal, 8 digits
        self.Total_salt = np.zeros(self.NumDays)              #Total mass of salt in the root zone, mg/L m, decimal, 8 digits      
        self.SaltDiff = np.zeros(self.NumDays)              #Difference in total salt between day before and present day (should equal sum of fluxes), mg/L m, decimal, 8 digits
        self.Ks_salt = np.zeros(self.NumDays)                #The reduction in ET due to salt stress, fraction, 5 digits
        self.Ks_nit = np.zeros(self.NumDays)                  #The reduction in ET due to nitrogen stress, fraction, 5 digits
        self.TotalIrrN = np.zeros(self.NumDays)         #Total nitrate added by irrigation to soil profile on each day, mg/L m, decimal, 8 digits
        self.TotalDrainN = np.zeros(self.NumDays)          #Total nitrate removed by drainage from soil profile on each day, mg/L m, decimal, 8 digits
        self.SumReactionsN = np.zeros(self.NumDays)       #Total of daily reactions (should equal sum of four rections), mg/L m, decimal, 8 digits
        self.SumFluxReactionsN = np.zeros(self.NumDays)  #Total of daily fluxes and reactions (should equal mass difference), mg/L m, decimal, 8 digits
        self.Total_Min = np.zeros(self.NumDays)          #Total mineralization in soil profile on each day, mg/L m, decimal, 8 digits
        self.Total_Den = np.zeros(self.NumDays)         #Total denitrification in soil profile on each day, mg/L m, decimal, 8 digits
        self.Total_Fer = np.zeros(self.NumDays)        #Total fertilization in soil profile on each day, mg/L m, decimal, 8 digits
        self.Total_Upt = np.zeros(self.NumDays)         #Total uptake in soil profile on each day, mg/L m, decimal, 8 digits
        self.UptakeReq = np.zeros(self.NumDays)         #Plant uptake requirement, mg/L m, decimal, 8 digit
        self.Total_nit = np.zeros(self.NumDays)          #Total nitrate in soil profile on each day, mg/L m, decimal, 8 digits
        self.Total_nit_accum = np.zeros(self.NumDays)          #Total nitrate accumulated in season, mg/L m, decimal, 8 digits
        self.N_ave = np.zeros(self.NumDays)                #Total mineralization in soil profile on each day, mg/L m, decimal, 8 digits
        self.N_optimal_low = np.zeros(self.NumDays)   #'Calculated with range_frac in N simulation section, lower limit of optimal soil N mg/kg
        self.N_optimal_high = np.zeros(self.NumDays) #'Calculated with range_frac in N simulation section, upper limit of optimal soil N mg/kg
        self.Cum_Min = np.zeros(self.NumDays)              #Cumulative mineralization in soil profile during season, mg/L m, decimal, 8 digits
        self.Cum_Den = np.zeros(self.NumDays)             #Cumulative denitrification in soil profile during season, mg/L m, decimal, 8 digits
        self.Cum_Fer = np.zeros(self.NumDays)             #Cumulative fertilization in soil profile during season, mg/L m, decimal, 8 digits
        self.Cum_Upt = np.zeros(self.NumDays)              #Cumulative uptake in soil profile during season, mg/L m, decimal, 8 digits
        self.CumIrrN = np.zeros(self.NumDays)             #Cumulative nitrate added by irrigation to soil profile during season, mg/L m, decimal, 8 digits
        self.CumDrnN = np.zeros(self.NumDays)             #Cumulative nitrate removed by drainage from soil profile during season, mg/L m, decimal, 8 digits
        self.CumChangeN = np.zeros(self.NumDays)       #Cumulative nitrate change in soil profile during season, mg/L m, decimal, 8 digits
        self.EC_leach_eqn = np.zeros(self.NumDays)
        self.Ks_water_total = np.zeros(self.NumDays)
        self.Kcmax = np.zeros(self.NumDays)
        self.Kr = np.zeros(self.NumDays)
        self.Few = np.zeros(self.NumDays)
        self.E_wet = np.zeros(self.NumDays)
        self.Rain_infilt = np.zeros(self.NumDays)
        self.dry_Rain_infilt = np.zeros(self.NumDays)
        self.wet_Rain_infilt = np.zeros(self.NumDays)
        self.ET_trans_wet = np.zeros(self.NumDays)
        self.Ps = np.zeros(self.NumDays)
        self.ECe_ave_effective = np.zeros(self.NumDays)
        self.Fert1 = np.zeros(self.NumDays)
        self.Fert2 = np.zeros(self.NumDays)
        self.Fert3 = np.zeros(self.NumDays)
        self.Fert = np.zeros(self.NumDays)
        self.Irr_Nit = np.zeros(self.NumDays)
        self.N_leach_eqn = np.zeros(self.NumDays)
        self.N_upt = np.zeros(self.NumDays)
        self.N_max = np.zeros(self.NumDays)
        self.mupt_entire_profile = np.zeros(self.NumDays)   
        self.Nitrogen_Kcb = np.zeros(self.NumDays)   
        self.KN_low = np.zeros(self.NumDays)   
        self.KN_high = np.zeros(self.NumDays)
        
        if self.SimStartDOY > self.StartDOY:
            self.DOY[0] = int(OutputArray['DOY'])
            self.NDVI[0] = OutputArray['NDVI']                        #The NDVI measured by remote sensing for the location, decimal 8 digits
            self.NDVI_regression[0] = OutputArray['NDVI_regression']  #The interpolated NDVI by regression, decimal 8 digits
            self.Kcb_NDVI[0] = OutputArray['Kcb_NDVI']                #The Kcb calculated from the NDVI regression, decimal 8 digits
            self.Kcb_calculated[0] = OutputArray['Kcb_calculated']    #The Kcb calculated from FAO equations, decimal 8 digits
            self.Kcb[0] = OutputArray['Kcb_used']                     #The Kcb used in the program, decimal 8 digits                            
            self.one_minus_c[0] = OutputArray['one_minus_c']          #The one minus the fraction of the canopy area (ground area fraction), decimal 8 digits
            self.Crop_height[0] = OutputArray['Crop_height']          #The h8 of the crop in meters, decimal 8 digits
            self.Root[0] = OutputArray['Root_depth']                  #The depth of the root in meters, decimal 8 digits
            self.Ke[0] = OutputArray['Ke']                            #The evaporation coefficient, decimal 8 digits
            self.E[0] = OutputArray['Evap']                           #The evaporation from the soil surface, mm, decimal 8 digits
            self.ET_trans[0] = OutputArray['Transpiration']           #The transpiration from the plant, mm, decimal 8 digits
            self.ET[0] = OutputArray['ET']                            #The evapotranspiration from the crop, mm, decimal 8 digits
            self.ET_pot[0] = OutputArray['Pottranspiration']          #The potential transpiration from the crop, decimal 8 digits
            self.ETcum[0] = OutputArray['CumulativePotentialET']      #The cumulative.tual evapotranspiration for the season, decimal 8 digits
            self.ET_refO[0] = OutputArray['ReferenceET']                 #The evapotranspiration from the weather station, decimal 8 digits
            self.Ky[0] = OutputArray['Ky']                            #The daily crop sensitivity to water stress
            self.Depletion_total[0] = OutputArray['Depletion_total']  #The depth of depletion, mm, decimal 8 digits
            self.Percent_depletion_total[0] = OutputArray['Percent_depletion_total'] #The percent depletion, percentage, decimal 8 digits
            self.Irrigation[0] = OutputArray['Irrigation']            #The depth of irrigation applied, mm, decimal 8 digits
            self.P[0] = OutputArray['P']                              #The management allowed depletion, fraction, decimal 8 digits
            self.Rainfall[0] = OutputArray['Rainfall']                #The rainfall depth, mm, decimal 8 digits
            self.Rain_Infilt[0] = OutputArray['Rain_infilt']
            self.FCave[0] = OutputArray['Field_capacity_ave']                #The average field capacity for soil profilel, decimal 8 digits
            self.PWPave[0] = OutputArray['Permanent_wilting_point_ave']      #The average permanent wilting point for soil profilel, decimal 8 digits
            self.TAWsum[0] = OutputArray['TAWsum']                          #Total available water, mm, decimal 8 digits
            self.RAWsum[0] = OutputArray['RAWsum']                          #Readily available water, mm, decimal 8 digits
            self.AW[0] = OutputArray['Available_water']        #The amount of water left in the soil for the plant, mm, decimal 8 digits
            self.VWCave[0] = OutputArray['Volumetric_water_content_ave']     #The average volumetric water content for the soil profile, fraction, decimal 8 digits
            self.TDW[0] = OutputArray['Total_depth_water']            #The total depth of water in the soil profile, mm, decimal 8 digits           
            self.zwt[0] = OutputArray['Water_table_elevation']         #The water table elevation above the datum, m, decimal 8 digits
            self.zu[0] = OutputArray['HeightEqu']                      #The upper elevation of the highest layer in equilibrium with the water table above the datum, m, 8 digits 
            self.Equilibrium_Max[0] = OutputArray['NumberEqu']         #The number of the upper layer that is in equilibrium with the water table, integer
            self.zu_Eq_Max = OutputArray['HeightEqu']
            self.Ks[0] = OutputArray['Water_stress']                #Ks which is the fraction of decreased evapotranspiration due to limited water, decimal, 8 digits
            self.Irr_Salt_Added[0] = OutputArray['IrrigationSaltAdded']      #The salt mass added by irrigation water, mg/L m, decimal, 8 digits
            self.Waste_Salt_Added[0] = OutputArray['WasteSaltAdded']         #The salt mass added by waste, mg/L m, decimal, 8 digits
            self.Seepage[0] = OutputArray['SeepageSaltLost']           #The salt mass lost by seepage, mg/L m, decimal, 8 digits
            self.SumOfFluxes[0] = OutputArray['SumOfFluxes']           #The mass balance of all salts added or lost, mg/L m, decimal, 8 digits
            self.Total_salt[0] = OutputArray['TotalSalt']               #Total mass of salt in the root zone, mg/L m, decimal, 8 digits      
            self.SaltDiff[0] = OutputArray['DifferenceInSalt']         #Difference in total salt between day before and present day (should equal sum of fluxes), mg/L m, decimal, 8 digits
            self.Ks_salt[0] = OutputArray['SaltStress']                #The reduction in ET due to salt stress, fraction, 5 digits
            self.Ks_nit[0] = OutputArray['NitStress']                  #The reduction in ET due to nitrogen stress, fraction, 5 digits
            self.TotalIrrN[0] = OutputArray['TotalIrrN']                 #Total nitrate added by irrigation to soil profile on each day, mg/L m, decimal, 8 digits
            self.TotalDrainN[0] = OutputArray['TotalDrainN']           #Total nitrate removed by drainage from soil profile on each day, mg/L m, decimal, 8 digits
            self.SumReactionsN[0] = OutputArray['SumReactionsN']       #Total of daily reactions (should equal sum of four rections), mg/L m, decimal, 8 digits
            self.SumFluxReactionsN[0] = OutputArray['SumFluxReactionsN']  #Total of daily fluxes and reactions (should equal mass difference), mg/L m, decimal, 8 digits
            self.Total_Min[0] = OutputArray['TotalMin']          #Total mineralization in soil profile on each day, mg/L m, decimal, 8 digits
            self.Total_Den[0] = OutputArray['TotalDen']          #Total denitrification in soil profile on each day, mg/L m, decimal, 8 digits
            self.Total_Fer[0] = OutputArray['TotalFer']          #Total fertilization in soil profile on each day, mg/L m, decimal, 8 digits
            self.Total_Upt[0] = OutputArray['TotalUpt']          #Total uptake in soil profile on each day, mg/L m, decimal, 8 digits
            self.UptakeReq[0] = OutputArray['UptakeReq']         #Plant uptake requirement, mg/L m, decimal, 8 digit
            self.Total_nit[0] = OutputArray['TotalNit']          #Total nitrate in soil profile on each day, mg/L m, decimal, 8 digits
            self.Total_nit_accum[0] = OutputArray['TotalNitAccum']          #Total nitrate in soil profile on each day, mg/L m, decimal, 8 digits
            self.N_ave[0] = OutputArray['NitAve']                #Total mineralization in soil profile on each day, mg/L m, decimal, 8 digits
            self.N_optimal_low[0] = OutputArray['NOptimalLow']   #'Calculated with range_frac in N simulation section, lower limit of optimal soil N mg/kg
            self.N_optimal_high[0] = OutputArray['NOptimalHigh'] #'Calculated with range_frac in N simulation section, upper limit of optimal soil N mg/kg
            self.Cum_Min[0] = OutputArray['CumMin']              #Cumulative mineralization in soil profile during season, mg/L m, decimal, 8 digits
            self.Cum_Den[0] = OutputArray['CumDen']              #Cumulative denitrification in soil profile during season, mg/L m, decimal, 8 digits
            self.Cum_Fer[0] = OutputArray['CumFer']              #Cumulative fertilization in soil profile during season, mg/L m, decimal, 8 digits
            self.Cum_Upt[0] = OutputArray['CumUpt']              #Cumulative uptake in soil profile during season, mg/L m, decimal, 8 digits
            self.CumIrrN[0] = OutputArray['CumIrrN']             #Cumulative nitrate added by irrigation to soil profile during season, mg/L m, decimal, 8 digits
            self.CumDrnN[0] = OutputArray['CumDrnN']             #Cumulative nitrate removed by drainage from soil profile during season, mg/L m, decimal, 8 digits
            self.CumChangeN[0] = OutputArray['CumChangeN']       #Cumulative nitrate change in soil profile during season, mg/L m, decimal, 8 digits
        self.DOY[0] = self.datelist[0].timetuple().tm_yday
        
        for j in range(1, self.NumDays - 1):
            if self.StartDOY >= self.SimStartDOY:
                adj_j = j + self.StartDOY
            else:
                adj_j = j + self.SimStartDOY
            if int(self.Select_Kcb) == 2:
                self.DOY[j] = self.datelist[j].timetuple().tm_yday
                self.KN_low[j]= 1.2
                self.KN_high[j] = 1.2
                if adj_j < 0:
                    self.Kcb[j] = 0
                elif adj_j < self.IP:
                    self.Kcb[j] = self.IKcb
                elif adj_j < self.IP + self.DP:
                    self.Kcb[j] = self.IKcb + (adj_j - self.IP) / (self.DP) * (self.MKcb - self.IKcb)
                elif adj_j < self.IP + self.DP + self.MP:
                    self.Kcb[j] = self.MKcb
                elif adj_j < self.IP + self.DP + self.MP + self.LP:
                    self.Kcb[j] = self.MKcb + (adj_j - (self.IP + self.DP + self.MP)) / (self.LP) * (self.EKcb - self.MKcb)
                elif adj_j < self.IP + self.DP + self.MP + self.LP + self.EP:
                    self.Kcb[j] = self.EKcb
                elif adj_j < self.IP + self.DP + self.MP + self.LP + self.EP + self.DP_2:
                    self.Kcb[j] = self.EKcb + (adj_j - (self.IP + self.DP + self.MP + self.LP + self.EP)) / (self.DP_2) * (self.MKcb_2 - self.EKcb)
                elif adj_j < self.IP + self.DP + self.MP + self.LP + self.EP + self.DP_2 + self.MP_2:
                    self.Kcb[j] = self.MKcb_2
                elif adj_j < self.IP + self.DP + self.MP + self.LP + self.EP + self.DP_2 + self.MP_2 + self.LP_2:
                    self.Kcb[j] = self.MKcb_2 + (adj_j - (self.IP + self.DP + self.MP + self.LP + self.EP + self.DP_2 + self.MP_2 + self.PD)) / (self.LP_2) * (self.EKcb_2 - self.MKcb_2)
                else:
                    self.Kcb[j] = self.EKcb_2
                    
#1 - c#1 - c#1 - c#1 - c#1 - c#1 - c#1 - c#1 - c#1 - c#1 - c#1 - c#1 - c
                if adj_j < 0:
                    self.one_minus_c[j] = 1
                elif adj_j < self.ICP:
                    self.one_minus_c[j] = self.IC
                elif adj_j < self.ICP + self.DCP:
                    self.one_minus_c[j] = self.IC - (adj_j - self.ICP) / self.DCP * (self.IC - self.MC)
                elif adj_j < self.ICP + self.DCP + self.MCP:
                    self.one_minus_c[j] = self.MC
                elif adj_j < self.ICP + self.DCP + self.MCP + self.ECP:
                    self.one_minus_c[j] = self.MC - (adj_j - self.ICP - self.DCP - self.MCP) / self.ECP * (self.MC - self.EC)
                else:
                    self.one_minus_c[j] = self.EC
#crop height
                if adj_j < 0:
                    self.Crop_height[j] = 0
                elif adj_j < self.ICP:
                    self.Crop_height[j] = self.minH
                elif adj_j < self.ICP + self.DCP:
                    self.Crop_height[j] = self.minH + (adj_j - self.ICP) / self.DCP * (self.midH - self.minH)
                elif adj_j < self.DCP + self.MCP + self.ICP:
                    self.Crop_height[j] = self.midH
                elif adj_j < self.DCP + self.MCP + self.ICP + self.ERH:
                    self.Crop_height[j] = self.midH + (adj_j - self.ICP - self.DCP - self.MCP) / self.ECP * (self.endH - self.midH)
                else:
                    self.Crop_height[j] = self.endH
#'roots
                if adj_j < 0:
                    self.Root[j] = 0
                elif adj_j < self.ICP:
                    self.Root[j] = self.IR * 1000
                elif adj_j < self.ICP + self.DCP:
                    self.Root[j] = (self.IR + (adj_j - self.ICP) / self.DCP * (self.DR - self.IR)) * 1000
                elif adj_j < self.DCP + self.MCP + self.ICP:
                    self.Root[j] = self.DR * 1000
                elif adj_j < self.DCP + self.MCP + self.ICP + self.ERH:
                    self.Root[j] = (self.DR + (adj_j - self.ICP - self.DCP - self.MCP) / self.ERH * (self.FR - self.DR)) * 1000
                else:
                    self.Root[j] = self.FR * 1000
                    
            self.Ps[1] = self.pTable22                
        
        if self.Salinity_simulation == 1:
            if self.WasteAppTF == 1:
                self.waste_app_day = (self.waste_app_date - self.DateP).days
            else:
                self.waste_app_day = 10000
            for j in range(1, self.NumDays):
                if self.StartDOY > self.SimStartDOY:
                    adj_j = j + self.StartDOY
                else:
                    adj_j = j + self.SimStartDOY
                self.Irr_Sal[j] = self.ECiw
                self.EC_leach_eqn[j] = 0
                if adj_j - self.waste_app_day > 0 and adj_j - self.waste_app_day - 1 < 1/self.Waste_dissolution:
                    self.Waste_sal[j] = self.Waste_dissolution * self.Waste_sal_kgha

        if self.Nitrogen_simulation == 1:
            if self.WasteAppTF == 1:
                self.Fert1_day = (self.Fert1_date - self.DateP).days
                self.Fert2_day = (self.Fert2_date - self.DateP).days
                self.Fert3_day = (self.Fert3_date - self.DateP).days
            else:
                self.Fert1_day = 10000
                self.Fert2_day = 10000
                self.Fert3_day = 10000

            self.Nitrogen_frac_sum = 0

            for j in range(1, self.NumDays):
                if self.StartDOY >= self.SimStartDOY:
                    adj_j = j + self.StartDOY
                else:
                    adj_j = j + self.SimStartDOY
                self.Fert[j] = 0
                self.Irr_Nit[j] = self.Niw
                self.N_leach_eqn[j] = 0
                if self.Number_fertilizations > 0 and adj_j - self.Fert1_day > 0 and adj_j - self.Fert1_day - 1 < 1/self.nit_dissolution:
                    self.Fert1[j] = self.nit_dissolution * self.Fert1_rate
                if self.Number_fertilizations > 1 and adj_j - self.Fert2_day > 0 and adj_j - self.Fert2_day - 1 < 1/self.nit_dissolution:
                    self.Fert2[j] = self.nit_dissolution * self.Fert2_rate
                if self.Number_fertilizations > 2 and adj_j - self.Fert3_day > 0 and adj_j - self.Fert3_day - 1 < 1/self.nit_dissolution:
                    self.Fert3[j] = self.nit_dissolution * self.Fert3_rate
                self.Fert[j] = self.Fert1[j] + self.Fert2[j] + self.Fert3[j]   
                
                # if adj_j < self.dev:
                if adj_j < self.IP:
                    self.Nitrogen_Kcb[j] = 0
                elif adj_j < self.IP + self.DP:
                    self.Nitrogen_Kcb[j] = (self.MKcb - self.IKcb) * ((adj_j - self.IP) / (self.IP + self.DP - self.IP))
                    # self.Nitrogen_Kcb[j] = (self.Kcb_mid - self.Kcb_initial) * ((adj_j - self.dev) / (self.mid - self.dev))
                elif adj_j < self.LP:
                    self.Nitrogen_Kcb[j] = self.MKcb - self.IKcb
                elif adj_j < self.EP:
                    self.Nitrogen_Kcb[j] = self.MKcb - self.IKcb - (self.MKcb - self.Kcb_end) * ((adj_j - self.EP) / (self.EP - self.LP))
                else:
                   self.Nitrogen_Kcb[j] = 0
                self.Nitrogen_frac_sum = self.Nitrogen_frac_sum + self.Nitrogen_Kcb[j]
                if self.Nitrogen_frac_sum == 0:
                    self.Nitrogen_frac_sum = 1

class output_layers:
    def __init__(self, Output_Layer_Array):
        self.DAP = np.zeros((self.NumDays, self.Num_layers + 2))                      #days after planting, integer
        self.PlantingIDL = np.zeros((self.NumDays, self.Num_layers + 2))              #The long random identifier that represents the planting, for example, 9a649ac3a4353ffe6d67fdad0c6bf3a9
        self.PlantingNameL = np.zeros((self.NumDays, self.Num_layers + 2))            #The name of the planting (text), normally includes owner name, field name, crop, and year, such as GaryField1Wheat2017
        self.LocationL = np.zeros((self.NumDays, self.Num_layers + 2))                #The position (location) in the planting, integer
        self.Layer = np.zeros((self.NumDays, self.Num_layers + 2))                    #The layer, integer        
        self.WC = np.zeros((self.NumDays, self.Num_layers + 2))                       #The water content in the soil layer, decimal, 8 digits
        self.Depletion = np.zeros((self.NumDays, self.Num_layers + 2))                #The depletion of water below field capacity (mm) in the soil layer, decimal, 8 digits
        self.Act_frac = np.zeros((self.NumDays, self.Num_layers + 2))                 #The fraction of total evaporated water that is removed from each soil layer, decimal, 8 digits
        self.Infilt = np.zeros((self.NumDays, self.Num_layers + 2))                   #The fraction of totol ET water that is evaporated from the particular layer, decimal, 8 digits
        self.Percent_depletion = np.zeros((self.NumDays, self.Num_layers + 2))        #The percent of the water between field capacity and permanent wilting point that is depleted from the layer, decimal, 8 digits
        self.EC = np.zeros((self.NumDays, self.Num_layers + 2))                       #The electrical conductivity of the actual water in the layer, dS/m, decimal, 8 digits
        self.ECe = np.zeros((self.NumDays, self.Num_layers + 2))                      #The electrical conductivity of the saturated paste extract from the layer, dS/m, decimal, 8 digits
        self.Mass_salt = np.zeros((self.NumDays, self.Num_layers + 2))                #The total mass of salt in the layer, decimal, 8 digits
        self.N = np.zeros((self.NumDays, self.Num_layers + 2))                        #Nitrate concentration in soil water, mg/L, decimal, 8 digits
        self.N_soil = np.zeros((self.NumDays, self.Num_layers + 2))                   #Nitrate concentration per mass of soil in each layer, mg/kg, decimal, 8 digits
        self.hc = np.zeros((self.NumDays, self.Num_layers + 2))                   #Nitrate concentration per mass of soil in each layer, mg/kg, decimal, 8 digits
        self.WC_test = np.zeros((self.Num_layers + 2, 102))
        self.WC_depth = np.zeros((self.NumDays, self.Num_layers + 2))
        self.Sum_depths = np.zeros(self.NumDays)
        self.Bottom_layer = np.zeros(self.NumDays)
        self.ts_depth = np.zeros(self.Num_layers + 2)
        self.Total_TAW_Fw = np.zeros((self.NumDays, self.Num_layers + 3))
        self.Uptake = np.zeros((self.NumDays, self.Num_layers + 2))
        self.fmntemp = np.zeros((self.NumDays, self.Num_layers + 2))
        self.kgha = np.zeros((self.NumDays, self.Num_layers + 2))
        self.Net_accumulation = np.zeros((self.NumDays, self.Num_layers + 2))
        self.Cum_Drain = np.zeros((self.NumDays, self.Num_layers + 2))
        self.Cum_Irr = np.zeros((self.NumDays, self.Num_layers + 2))
        self.Cum_Total = np.zeros((self.NumDays, self.Num_layers + 2))
        
        if len(Output_Layer_Array[Output_Layer_Array['Layer']==1].index.values) > 0:
            kk = int(Output_Layer_Array[Output_Layer_Array['Layer']==1].index.values)
        
        for k in range(1, self.Num_layers + 2):
            if self.SimStartDOY <= self.StartDOY:   
                self.WC[0][k] = self.InitWC[k]
                print('initial water content',self.InitWC[k])
                self.Depletion[0][k] = (self.FC[k] - self.WC[0][k]) * self.dz[k] * self.FW_layers[1][k]
                self.Percent_depletion[0][k] = (self.FC[k] - self.WC[0][k]) / (self.FC[k] - self.PWP[k])
            else:
                self.WC[0][k] = Output_Layer_Array['WaterContent'].loc[kk]
                self.Depletion[0][k] = Output_Layer_Array['Depletion'].loc[kk]
                self.Percent_depletion[0][k] = Output_Layer_Array['PercentDepletion'].loc[kk]
        for j in range(1, self.NumDays - 1):
            self.Total_TAW_Fw[j][self.Num_layers + 2] = 0
            for k in range(self.Num_layers + 1, 0, -1):
                self.Total_TAW_Fw[j][k] = self.Total_TAW_Fw[j][k+1] + self.TAW[k]*self.FW_layers[j][k]            
            k = self.Num_layers + 1
            while self.Depth[k] < self.Root[j]/1000:
                k = k - 1
            self.Bottom_layer[j] = k
            if self.Bottom_layer[j] < 1:
                self.Bottom_layer[j] = 1
            self.FCave[j] = 0
            self.PWPave[j] = 0
            self.TAWsum[j] = 0
            for k in range(int(self.Bottom_layer[j]), self.Num_layers + 2):
                self.FCave[j] = self.FC[k] + self.FCave[j]    
                self.PWPave[j] = self.PWP[k] + self.PWPave[j]  
                self.TAWsum[j]= self.TAWsum[j] + (self.FC[k] - self.PWP[k]) * self.dz[k]
            self.RAWsum[j] = self.TAWsum[j] * self.Ps[j]
            self.FCave[j] = self.FCave[j] / (self.Num_layers - self.Bottom_layer[j] + 1)
            self.PWPave[j] = self.PWPave[j] / (self.Num_layers - self.Bottom_layer[j] + 1)
        sum_ET = 0
        for k in range(int(self.Bottom_layer[1]), self.Num_layers + 2):
            sum_ET = sum_ET + self.ET_fractions[1][k]

        if self.No_ET_frac_Adjustment == 1:
            for k in range(1, self.Num_layers + 2):
                self.Act_frac[1][k] = self.ET_fractions[1][k] / sum_ET
        else:
            Remaining = 1
            Sum_act_frac = 0
            for k in range(self.Num_layers + 1, int(self.Bottom_layer[1]) - 1, -1):
                self.RAW[k] = self.Ps[1] * self.TAW[k]
                if self.Fw[k] < self.Fw_y[k]:
                    self.Depletion[0][k] = self.Depletion[0][k] * self.Fw[k] / self.Fw_y[k]
                if (self.Depletion[0][k] > self.RAW[k] * self.Fw[k]) and (self.Depletion[0][k] < self.TAW[k] * self.Fw[k]):
                    self.Act_frac[1][k] = self.ET_fractions[j][k] * (self.TAW[k] * self.Fw[k] - self.Depletion[0][k]) / ((self.TAW[k] - self.RAW[k]) * self.Fw[k])/sum_ET
                elif self.Depletion[0][k] > self.TAW[k] * self.Fw[k]:
                    self.Act_frac[1][k] = 0
                else:
                    self.Act_frac[1][k] = self.ET_fractions[1][k]/sum_ET
                Sum_act_frac = Sum_act_frac + self.Act_frac[1][k]
            
            if Sum_act_frac > 0:
                for k in range(self.Num_layers + 1, int(self.Bottom_layer[1]) - 1, -1):
                    self.Act_frac[1][k] = self.Act_frac[1][k] / Sum_act_frac

        if self.Salinity_simulation == 1:
            for k in range(1, self.Num_layers + 2):
                kk = Output_Layer_Array[Output_Layer_Array['Layer']==k].index.values
                if self.SimStartDOY <= self.StartDOY:                                                          #The day of experiment, integer j
                    self.ECe[0][k] = self.InitECe[k]
                    self.EC[0][k] = self.ECe[0][k] / self.WC[0][k] * self.ts[k]  
                else:
                    self.ECe[0][k] = Output_Layer_Array['ECe'].loc[kk]
                    self.EC[0][k] = Output_Layer_Array['EC'].loc[kk] 
            Act_frac_sum = 0
            for k in range(self.Num_layers, int(self.Bottom_layer[1]) - 1, -1):
                self.ECe_ave_effective[1] = self.ECe_ave_effective[1] + self.Act_frac[1][k] * self.ECe[0][k] * self.FC[k] / self.WC[0][k]
                Act_frac_sum = Act_frac_sum + self.Act_frac[1][k]
            if Act_frac_sum == 0:
                Act_frac_sum = 1
            self.ECe_ave_effective[1] = self.ECe_ave_effective[1] / Act_frac_sum
            self.Total_salt[0] = 0
            for k in range(1, self.Num_layers + 2):
                self.Total_salt[0] = self.Total_salt[0] + self.EC[0][k] * 640 * self.WC[0][k] * self.dz[k] * self.Fw[k]
            for k in range(1, self.Num_layers + 2):
                self.Active_EC[k] = self.EC[0][k]
                                    #The electrical conductivity of the actual water in the layer, dS/m, decimal, 8 digits
        if self.Nitrogen_simulation == 1:
            for k in range(1,  self.Num_layers + 2):
                kk = Output_Layer_Array[Output_Layer_Array['Layer']==k].index.values
                if self.SimStartDOY <= self.StartDOY:                                                          #The day of experiment, integer j
                    self.N_soil[0][k] = self.InitsoilN[k]       #Nitrate concentration per mass of soil in each layer, mg/kg, decimal, 8 digits
                    self.N[0][k] = self.N_soil[0][k] * self.bd[k] /self.WC[0][k]                #Nitrate concentration in soil water, mg/L, decimal, 8 digits
                    self.kgha[0][k] = self.N[0][k] * self.dz[k] * 10 * self.WC[0][k]
                else:
                    self.N_soil[0][k] = Output_Layer_Array['Nitratemgkg'].loc[kk]       #Nitrate concentration per mass of soil in each layer, mg/kg, decimal, 8 digits
                    self.N[0][k] = Output_Layer_Array['Nitrate'].loc[kk]               #Nitrate concentration in soil water, mg/L, decimal, 8 digits
                    self.kgha[0][k] = self.N[0][k] * self.dz[k] * 10 * self.WC[0][k]
            for i in range(1, self.NumDays):
                self.N_upt[i] = self.Seasonal_N_uptake * self.Nitrogen_Kcb[i] / self.Nitrogen_frac_sum
                temp = self.T_bar + self.A0 * np.sin(2 * 3.14159 / 365 * (self.DOY[i] - self.tm))
                for k in range(1, self.Num_layers + 1):
                    self.fmntemp[i][k] = self.Qtemp ** ((temp - 20) / 10)
                    if self.fmntemp[i][k] > 1.5:
                        self.fmntemp[i][k] = 1.5 - (self.fmntemp[i][k] - 1.5) / 5

    def Create_output_layer_arrays(self, Num, Num_days, Num_layers):
        self.date_layer_out = np.empty(Num, dtype = dt)
        self.Days_after_planting_layer = np.zeros(Num)
        self.pid_layer_out  = np.empty(Num, dtype = object)
        self.planting_name_layer_out =  np.empty(Num, dtype = object)
        self.location_layer_out = np.zeros(Num)
        self.layer_out = np.zeros(Num)
        self.WC_layer_out = np.zeros(Num)
        self.EC_out = np.zeros(Num)
        self.N_out = np.zeros(Num)
        self.Depletion_layer_out = np.zeros(Num)
        self.ActFrac_layer_out = np.zeros(Num)
        self.Infiltration_layer_out =  np.zeros(Num)
        self.Percent_depletion_layer_out = np.zeros(Num)
        self.EC_layer_out = np.zeros(Num)
        self.ECe_layer_out = np.zeros(Num)
        self.Mass_salt_layer_out = np.zeros(Num)
        self.N_layer_out = np.zeros(Num)
        self.N_soil_layer_out = np.zeros(Num)
        
        l = 0
        for j in range(1, Num_days):
            for k in range(1, Num_layers + 2):
                l = l + 1
                self.date_layer_out[l] = self.datelist[j]
                self.Days_after_planting_layer[l] = self.DOE[j]
                self.pid_layer_out[l]  = self.PlantingIDO[j]
                self.planting_name_layer_out[l] = self.PlantingNameO[j]
                self.location_layer_out[l] = self.LocationO[j]
                self.layer_out[l] = k
                self.WC_layer_out[l] = self.WC[j][k]
                self.Depletion_layer_out[l] = self.Depletion[j][k]
                self.ActFrac_layer_out[l] = self.Act_frac[j][k]
                self.Infiltration_layer_out[l] =  self.Infilt[j][k]
                self.Percent_depletion_layer_out[l]
                self.EC_layer_out[l] = self.EC[j][k]
                self.ECe_layer_out[l] = self.ECe[j][k]
                self.Mass_salt_layer_out[l] = self.Mass_salt[j][k]
                self.N_layer_out[l] = self.N[j][k]
                self.N_soil_layer_out[l] = self.N_soil[j][k]



class model(plantings, weather, fields, status, soil, ET_fractions, Wetting_fractions, output, output_layers, irrigation, ET_daily, RS_daily):
    def __init__(self, Planting_Array,
                 Weather_Array, 
                 Field_Array, 
                 Status_Array,
                 Soil_Array,
                 ETfrac_Array, 
                 Wetting_Array, 
                 Output_Array,
                 Output_Layers_Array,
                 Irrigation_Array,
                 ET_Daily_Array,
                 RS_Daily_Array):
        plantings.__init__(self, Planting_Array)
        weather.__init__(self, Weather_Array)
        fields.__init__(self, Field_Array)
        status.__init__(self, Status_Array)
        soil.__init__(self, Soil_Array)
        ET_fractions.__init__(self, ETfrac_Array)
        Wetting_fractions.__init__(self, Wetting_Array)
        output.__init__(self,Output_Array)
        output_layers.__init__(self, Output_Layers_Array)
        irrigation.__init__(self, Irrigation_Array)
        ET_daily.__init__(self, ET_Daily_Array)
        RS_daily.__init__(self, RS_Daily_Array)
        
       
    def wetting_front(self, j):
        m = np.zeros(self.Num_layers + 2) 
        Depth = np.zeros(self.Num_layers + 2)
        k = np.zeros(self.Num_layers + 2)
        z = np.zeros(self.Num_layers + 2)
        SAV = np.zeros(self.Num_layers + 2)
        Precip_rate = self.rainfall[j] / self.Rain_time[j] / 1000
        for i in range(1, self.Num_layers + 1):
            SAV[i] = self.SAV
            m[i] = self.ts[self.Num_layers + 1 - i] - self.WC[j][self.Num_layers + 1 - i]
            k[i] = self.Ksat[self.Num_layers + 1 - i] / 24
            z[i] = self.dz[self.Num_layers + 1 - i]
        
        Depth[0] = 0
        Capacity = 0
        Lf_prev = 0
        
        for i in range(1, self.Num_layers + 1):
            Depth[i] = Depth[i - 1] + z[i]
            Capacity = Capacity + z[i] * m[i]
        
        i = 1
        if m[i] < 0.01:
            while m[i] < 0.02 and i < self.Num_layers:
                i = i + 1
                Lf_prev = Depth[i]
        
        if Capacity > 0 and i < self.Num_layers:
            continuecalcs = 1
        else:
            continuecalcs = 0
            Total_Infiltration = Capacity
        
        total_time = 0
        Total_Infiltration = 0
        dt = self.Rain_time[j] / 100
        if Lf_prev > 0:
            LF = Lf_prev
        elif self.PondedDepth > 0:
            LF = 0.02
            Total_Infiltration = LF * m[1]
        else:
            LF = Precip_rate * dt / m[1]
            Total_Infiltration = LF * m[1]
        
        i = 1
        while Depth[i] < LF:
            i = i + 1
    
        j = i
        total_time = dt
        
        while total_time < self.Rain_time[j] and continuecalcs == 1 and j < self.Num_layers:
            total_time = total_time + dt
            Sum_denominator = 0
            for i in range(1, j+1):
                if i == j:
                    Sum_denominator = Sum_denominator + (LF - Depth[i - 1]) / k[i]
                else:
                    Sum_denominator = Sum_denominator + z[i] / k[i]
            Keff = LF / Sum_denominator
            Infiltration_rate = Keff * ((self.PondedDepth + self.SAV + LF) / LF)
            if Precip_rate < Infiltration_rate and self.PondedDepth > 0:
                Infiltration = Precip_rate * dt
            else:
                Infiltration = Infiltration_rate * dt
    
            Total_Infiltration = Total_Infiltration + Infiltration
            while (m[j] == 0 & continuecalcs == 1):
                j = j + 1
                LF = LF + z[j - 1]
                if j == self.Num_layers:
                    continuecalcs = 0
                    Total_Infiltration = Capacity
            if continuecalcs == 1:
                LF = LF + Infiltration / m[j]
                if LF > Depth[self.Num_layers]:
                    continuecalcs = 0
                    Total_Infiltration = Capacity
                else:
                    i = 1
                    while Depth[i] < LF:
                        i = i + 1
                    j = i
                    
            if self.Rain_Infilt[j] > self.rainfall[j]:
                self.Rain_Infilt[j] = self.rainfall[j]
            
        return Total_Infiltration;
        
    def parameterize_drainage(self, j):
        if j == 1:
            self.zwt[j - 1] = self.Init_zwt
            self.Equilibrium_Max[j - 1] = self.Equil_max_init
        
        self.ts_depth[0] = 0
        for i in range(1, self.Num_layers + 1):
            ts_depth[i] = ts_depth(i - 1) + self.ts[i] * self.dz[i]
            for ij in range(0, 101):
                self.zwt_test[i][ij] = ij * (self.zu[i]) / 100
                self.WC_test[i][ij] = 0
                for k in range(1, i + 1):
                    if self.zwt_test[i][ij] > self.zu[k]: #Water table is above cell
                        self.WC_test[i][ij] = self.ts[k] * self.dz[k] + self.WC_test[i][ij]
                    elif self.zu[k] > self.zwt_test[i][ij] and self.zl[k] < self.zwt_test[i][ij]:   #Water table is within cell
                        self.WC_test[i][ij] = d_cell_WT(self.ts[k], self.tr[k], self.zu[k], self.zl[k], self.zwt_test[i][ij], self.av[k], nv[k]) + self.WC_test[i][ij]
                    elif self.zu[k] == self.zwt_test[i][ij]: #Water table is at top of cell
                        self.WC_test[i][ij] = self.ts[k] * self.dz[k] + self.WC_test[i][ij]
                    else:    #Water table is below cell or at lower boundary of cell
                        self.WC_test[i][ij] = d_cell_Eq(self.ts[k], self.tr[k], self.zu[k], self.zl[k], self.zwt_test[i][ij], self.av[k], nv[k]) + self.WC_test[i][ij]
        
        
    def find_water_table_depths(self, j):    
        self.Sum_WC_depth[j - 1][0] = 0
        
        for k in range(1, self.Num_layers + 2):
            if zwt[j - 1] > self.zu[k]:
                self.WC[j - 1][k] = self.ts[k]
            elif self.zu[k] > zwt[j - 1] and self.zl[k] < zwt[j - 1]:
                self.WC[j - 1][k] = d_cell_WT(self.ts[k], self.tr[k], self.zu[k], self.zl[k], zwt[j - 1], self.av[k], nv[k]) / self.dz[k]
            elif self.zu[k] == zwt[j - 1]:
                self.WC[j - 1][k] = self.ts[k]
            else:
                self.WC[j - 1][k] = d_cell_Eq(self.ts[k], self.tr[k], self.zu[k], self.zl[k], zwt[j - 1], self.av[k], nv[k]) / self.dz[k]
    
            self.Sum_WC_depth[j - 1][k] = self.Sum_WC_depth[j - 1][k - 1] + self.WC[j - 1][k] * self.dz[k]
            Depletion[j - 1][k] = (self.FC[k] - self.WC[j - 1][k]) * self.dz[k] * self.Fw[k]
        self.WC[j - 1][0] = self.ts[1]
        self.Equilibrium_Max[j - 1] = 0
        self.Sum_WC_depth[j - 1][0] = 0
        for k in range(1, self.Num_layers + 2):
            self.Sum_WC_depth[j - 1][k] = self.Sum_WC_depth[j - 1][k - 1] + self.WC[j - 1][k] * self.dz[k]
            if self.Sum_WC_depth[j - 1][k] > ts_depth[k] * Fraction_of_saturation_for_equilibrium:
                self.Equilibrium_Max[j - 1] = k
    
        Sum_depths[j - 1] = 0
        if j > 1:
            Sum_depths[j - 2] = 0
        for k in range(1, self.Num_layers + 2):
            Sum_depths[j - 1] = Sum_depths[j - 1] + self.WC[j - 1][k] * self.dz[k] * self.Fw[k]
            if j > 1:
                Sum_depths[j - 2] = Sum_depths[j - 2] + self.WC[j - 1][k] * self.dz[k] * self.Fw[k]
    
    def water_table_simulation_with_eq_Max_greater_zero(self, j):
        ET_losses = 0
        for k in range(1, self.Equilibrium_Max[j] + 1):
            if k == self.Num_layers + 1:
                ET_losses = ET_losses + (self.E_wet[j] + self.ET_trans_wet[j] * self.Act_frac[j][k]) / 1000
            else:
                ET_losses = ET_losses + (self.ET_trans_wet[j] * self.Act_frac[j][k]) / 1000
        
#Find the drainage rate based on either the linear function or the Kirkham resistance function
        
        time_step = 1
        if self.Kirkham_rate == 1:
            if self.zwt[j-1] == self.Drain_elevation:
                Drain_flux = 0
            else:
                Drain_flux = self.Keff_horizontal / (1 + self.L_drain * self.Kirkham_F / (self.zwt[j-1] - self.Drain_elevation))
        else:
            if self.zwt[j-1] == self.Drain_elevation:
                Drain_flux = 0
            else:
                Drain_flux = self.Drain_rate * (self.zwt[j-1] - self.Drain_elevation) * time_step
        if Drain_flux < 0:
            Drain_flux = 0
#Find the fluxes above and below the range in equilibrium with the water table
        
        if self.Equilibrium_Max[j] == self.Num_layers + 1:
            Fluxes = self.Irrigation[j] / 1000 + self.Rain_Infilt[j] / 1000 - Drain_flux
        else:
            Fluxes = self.Infilt[j][self.Equilibrium_Max[j] + 1] - Drain_flux  #took out + 1
        
        
        k = self.Equilibrium_Max[j]
        
#Add the fluxes and ET losses to the sum of water depth from previous time step
#in the layers in equilibrium with the water table

        self.Sum_WC_depth[j][k] = self.Sum_WC_depth[j-1][k] + (Fluxes - ET_losses)
        
#Find the elevation of the water table that would result in the Sum_WC_depth depth of water in the
#equilibrium layers
                   
        i = 1   #i is the counter for WC_test, there are 100 possible elevations
                #i = 1 corresponds to a water table at the bottom of the soil profile
        while i < 101 and self.WC_test[k][i] < self.Sum_WC_depth[j][k]: #test begins in the
            i = i + 1
        
#if the depth is not found in the current layer corresponding to equilibrium_max
#then the program looks in the next layer up and incorporates a layer above the equilibrium_max level
        
        Ave_WC = self.Sum_WC_depth[j][k] / (zu(self.Equilibrium_Max[j]) - self.zl[1])
        
#Set the water table at the elevation found, otherwise set the water table at the top of the soil profile
        
        if i < 101:
            self.zwt[j] = zwt_test(k, i)
        else:
            self.zwt[j] = self.zu(self.Num_layers + 1)
        
        
        self.zwt[j] = Find_wt(self.ts[k], self.tr[k], self.zu(self.Equilibrium_Max[j]), self.zl[1], self.av[k], self.nv[k], Ave_WC, self.zwt[j])
        
        for k in range(1, self.Equilibrium_Max[j] + 1):
            if self.zwt[j] > self.zu[k]:  #Water table is above cell
                self.WC[j][k] = self.ts[k]
            elif self.zu[k] > self.zwt[j] and self.zl[k] < self.zwt[j]:   #Water table is within cell
                self.WC[j][k] = d_cell_WT(self.ts[k], tr[k], self.zu[k], self.zl[k], self.zwt[j], self.av[k], self.nv[k]) / self.dz[k]
            elif self.zu[k] == self.zwt[j]:  #Water table is at top of cell
                self.WC[j][k] = self.ts[k]
            else:    #Water table is below cell or at lower boundary of cell
                self.WC[j][k] = d_cell_Eq(self.ts[k], tr[k], self.zu[k], self.zl[k], self.zwt[j], self.av[k], self.nv[k]) / self.dz[k]
            
            
            if k == 1:
                self.Infilt[j][k] = Drain_flux
                self.WC[j][0] = self.WC[j][1]
            else:
                self.Infilt[j][k] = self.Infilt[j][k - 1] + ET_trans[j] / 1000 * self.Act_frac[j][k - 1] + (self.WC[j][k - 1] - self.WC[j - 1][k - 1]) * self.dz(k - 1)
            
            Depletion[j][k] = (self.FC[k] - self.WC[j][k]) * self.dz[k] * self.Fw[k]
            Percent_depletion[j][k] = (self.FC[k] - self.WC[j][k]) / (self.FC[k] - PWP[k]) * 100
        
        base = 0
        self.Equilibrium_Max[j] = 0
        self.Sum_WC_depth[j][0] = 0
        if Keep_equilibrium_layers_out_of_root_zone:
            Max_layer = Max_layer_of_equilibrium_layers
        else:
            Max_layer = self.Num_layers + 1
        

        for k in range(1, Max_layer + 1):
            self.Sum_WC_depth[j][k] = self.Sum_WC_depth[j][k - 1] + self.WC[j][k] * self.dz[k]
            if self.Sum_WC_depth[j][k] > ts_depth[k] * Fraction_of_saturation_for_equilibrium:
                self.Equilibrium_Max[j] = k
    def water_table_simulation_with_eq_Max_equals_zero(self, j):
        if Continue_drainage_rate:
            WC_stop = d_cell_Eq(self.ts[1], self.tr[1], self.zu[1], self.zl[1], 0, self.av[1], self.nv[1]) / self.dz[1]
            if self.WC[j - 1][1] < WC_stop:
                self.Infilt[j][1] = 0
            else:
                self.Infilt[j][1] = self.Infilt[j - 1][1] * self.Infilt[j - 1][1] / self.Infilt(j - 2, 1)
            
        else:
            self.Infilt[j][1] = 0
        
        self.WC[j][1] = self.WC[j - 1][1] - self.Infilt[j][1] / self.dz[1] / self.Fw[1] + self.Infilt[j][2] / self.dz[2] / self.Fw[2]
                
        
def hg(a, b, C, z):

    value = 1
    hg = value
    
    for i in range(1,41):
        value = value * (a + i - 1) * (b + i - 1) / (C + i - 1) / i * z
        hg = hg + value
    return hg

def d_total(ts, tr, zt, zwt, alpha, N):

    W = -((alpha * (zt - zwt) * 100) ** N)
    zwt = zwt
    hyper = hg(1, 1 - 1 / N, 1 + 1 / N, W / (W - 1))

    d_total = ts * zwt + tr * (zt - zwt) + (ts - tr) * (zt - zwt) * (1 - W) ** (1 / N - 1) * hyper
    return d_total

def d_cell_WT(ts, tr, zu, zl, zwt, alpha, N):

    W = -((alpha * (zu - zwt) * 100) ^ N)
    hyper = hg(1, 1 - 1 / N, 1 + 1 / N, W / (W - 1))

    d_cell_WT = ts * (zwt - zl) + tr * (zu - zwt) + (ts - tr) * (zu - zwt) * (1 - W) ** (1 / N - 1) * hyper
    return d_cell_WT

def d_cell_Eq(ts, tr, zu, zl, zwt, alpha, N):

    WU = -((alpha * (zu - zwt) * 100) ** N)
    wl = -((alpha * (zl - zwt) * 100) ** N)
    hyperu = hg(1, 1 - 1 / N, 1 + 1 / N, WU / (WU - 1))
    hyperl = hg(1, 1 - 1 / N, 1 + 1 / N, wl / (wl - 1))

    d_cell_Eq = tr * (zu - zl) + (ts - tr) * (zu - zwt) * (1 - WU) ** (1 / N - 1) * hyperu - (ts - tr) * (zl - zwt) * (1 - wl) ** (1 / N - 1) * hyperl
    return d_cell
    
def Find_wt(ts, tr, zu, zl, alpha, N, th_ave, zwt_last_time):

    dz = zu - zl
    d_layer_ave = th_ave * dz
    zwt = zwt_last_time
    
    W = -((alpha * (zu - zwt) * 100) ** N)
    
    hyper = hg(1, 1 - 1 / N, 1 + 1 / N, W / (W - 1))

    d_layer_2 = ts * (zwt - zl) + tr * (zu - zwt) + (ts - tr) * (zu - zwt) * (1 - W) ** (1 / N - 1) * hyper
    
    if zwt > 0.15:
        Convergence = 0.00001
    else:
        Convergence = 0.00001
   
    while abs(d_layer_2 - d_layer_ave) > Convergence:
    
        W = -((alpha * (zu - zwt) * 100) ** N)
        
        hyper = hg(1, 1 - 1 / N, 1 + 1 / N, W / (W - 1))
    
        d_layer_2 = ts * (zwt - zl) + tr * (zu - zwt) + (ts - tr) * (zu - zwt) * (1 - W) ** (1 / N - 1) * hyper
    
        zwt = zwt * (d_layer_ave / d_layer_2) ^ 0.5
    
    return zwt
        
def zwt_crit(ts, tr, zu, zl, zwt, alpha, N):

    W = -((alpha * (zu - zwt) * 100) ** N)
    hyper = hg(1, 1 - 1 / N, 1 + 1 / N, W / (W - 1))

    zwt_crit = ts * (zwt - zl) + tr * (zu - zwt) + (ts - tr) * (zu - zwt) * (1 - W) ** (1 / N - 1) * hyper
    return zwt_crit

def Calculate_k(te, mv, Ko, L):
        
    if te > 1:
        Calculate_k = Ko * 0.1
    else:
        Calculate_k = Ko * te ** L * (1 - (1 - te ** (1 / mv)) ** mv) ** 2
    return Calculate_k

def Calculate_te(WC, tr, ts):

    if WC > tr:
        Calculate_te = (WC - tr) / (ts - tr)
    else:
        Calculate_te = 0.00001
    return Calculate_te

def KirkhamSolution(x, S, r, h, k, re):
    
    Kirkham_height = (math.log(math.sin(3.14159 * x / (2 * S)) / math.sin(Pi * r / (2 * S))))
    sum = 0
    for m in range(1, 400):
        sum = sum + ((1 / m) * (math.cos(m * 3.14159 * r / S) - math.cos(m * 3.14159 * x / S)) * (math.exp(-m * 3.14159 * h / S) / np.sinh(m * 3.14159 * h / S)))
    Kirkham_height = Kirkham_height + sum

    Kirkham_height = re * 2 * S * Kirkham_height / (k * 3.14159)
    return Kirkham_height  



# def Percent_depletion():
#     self.Planting= np.zeros(100)
#     self.start_day
#     self.slope = np.zeros(10)
#     self.yint = np.zeros(10)
#     self.Runit = np.zeros(100)
#     self.Leap_year As Boolean
#     self.Days_to_irrigate = np.zeros(10)
#     self.Irrigation_name As String
#     Leap_year = Worksheets("Main").Range("C12").value
#     Number_of_sheets = Worksheets("Spatial_data").Range("K5").value
#     Num_treatments = Worksheets("Spatial_data").Range("K6").value
#     Num_Cells = Worksheets("Spatial_data").Range("K7").value
#     Num_plots = Worksheets("Spatial_data").Range("K26").value
#     self.data As Variant, days As Variant
#     self.Perform_analysis = np.zeros(100) As Boolean
#     self.Perform_plot_analysis = np.zeros(144) As Boolean
#     self.Depletions = np.zeros(100) As Single
#     Num_Sections = Worksheets("Spatial_data").Range("K6").value
#     self.Sort_depletions = np.zeros(16) As Single
#     self.Low_sequence = np.zeros(16) As Single
#     self.Average_depletion = np.zeros(1200) As Single
#     self.Total_in_quartile = np.zeros(16) As Single
#     self.Average_in_quartile = np.zeros(16, 1200) As Single
#     self.Average_mm_in_quartile = np.zeros(16, 1200) As Single
#     self.Total_mm_in_quartile = np.zeros(16) As Single
#     self.Ave_irr_required = np.zeros(1200) As Single
#
#     self.PMAD = np.zeros(10, 16, 1200) As Single
#     self.mm_depletion = np.zeros(10, 16, 1200) As Single
#     Planting_date_in_spatial_data = Worksheets("Spatial_data").Range("P9").value
#
#     For i = 1 To Num_treatments
#         Perform_analysis(i) = Worksheets("Spatial_data").Cells(i + 10, 13).value
#         If Perform_analysis(i) = 1 Then
#             Irrigation_name = Worksheets("Spatial_data").Cells(1 + i, 3).Formula
#             Irrigation_number = Right(Irrigation_name, 2)
#         End If
#     Next i
#
#     Crop_data_planting_days = Worksheets("Crop_data").Range("A3:IZ3")
#
#     For i = 1 To Number_of_sheets
#         If Planting_date_in_spatial_data Then
#             Planting_day = Worksheets("Spatial_data").Cells(10 + CInt((i - 8) / 16) + 1, 16).value
#         Else
#             Planting_day = Crop_data_planting_days(1, i + 1)
#         End If
#         Current_day_of_year = Worksheets("Main").Range("E3").value
#         Number_of_rows = Current_day_of_year - Planting_day + 2
#         Cell_range_percent_depletion = "IR3:IR" + CStr(Number_of_rows + 3)
#         Cell_range_mm_depletion = "IQ3:IQ" + CStr(Number_of_rows + 3)
#         Cell_range_days = "W3:W" + CStr(Number_of_rows + 3)
#         Cell_range_days_right = "W5:W" + CStr(Number_of_rows + 3)
#         Cell_range_pmad = "IT3:IT" + CStr(Number_of_rows + 3)
#         Cell_range_days_out = "A5:A" + CStr(Number_of_rows + 4)
#         Cell_range_right_days_out = "AB7:AB" + CStr(Number_of_rows + 4)
#         Cell_range_date_out = "AC7:AC" + CStr(Number_of_rows + 4)
#         Plot_number = Worksheets("Spatial_data").Cells(1 + i, 6).value
#         Irrigation_name = Worksheets("Spatial_data").Cells(1 + i, 3).Formula
#         Irrigation_number = Right(Irrigation_name, 2)
#         Irrigation_depths = Worksheets("Irr_" + Irrigation_number).Range("B1:B400")
#         Irrigation_summary_sheet = CSng(Irrigation_number)
#         Section_number = Worksheets("Spatial_data").Cells(1 + i, 5).value
#
#         If Perform_analysis(Irrigation_summary_sheet) = 1 Then
#             Worksheets("Irrigation_summary_" + Irrigation_number).Range("A5:AC1000").ClearContents
#             Worksheet_name = "C_" + CStr(Plot_number)
#             days = Worksheets(Worksheet_name).Range(Cell_range_days)
#             days_right = Worksheets(Worksheet_name).Range(Cell_range_days_right)
#             Percent_depletion_transfer = Worksheets(Worksheet_name).Range(Cell_range_percent_depletion)
#             PMAD_transfer = Worksheets(Worksheet_name).Range(Cell_range_pmad)
#             mm_depletion_transfer = Worksheets(Worksheet_name).Range(Cell_range_mm_depletion)
#             Worksheets("Irrigation_summary_" + Irrigation_number).Range(Cell_range_days_out) = days
#             Worksheets("Irrigation_summary_" + Irrigation_number).Range(Cell_range_right_days_out) = days_right
#             First_cell_for_date = Worksheets("Irrigation_summary_" + Irrigation_number).Range("AB7")
#             Efficiency = Worksheets("Irrigation_summary_" + Irrigation_number).Range("AF7")
#             For ii = 1 To Number_of_rows + 10
#                 Place_date = Worksheets("DOY dates").Cells(First_cell_for_date + ii + 1, 2).value
#                 Worksheets("Irrigation_summary_" + Irrigation_number).Cells(ii + 6, 29).value = Place_date
#                 Worksheets("Irrigation_summary_" + Irrigation_number).Cells(ii + 6, 30).value = Irrigation_depths(First_cell_for_date + ii, 1)
#             Next ii
#
#             For j = 1 To Number_of_rows + 1
#                 Worksheets("Irrigation_summary_" + Irrigation_number).Cells(j + 4, Section_number + 1).value = Percent_depletion_transfer(j, 1)
#                 mm_depletion(Irrigation_summary_sheet, Section_number, j) = mm_depletion_transfer(j, 1)
#                 PMAD(Irrigation_summary_sheet, Section_number, j) = PMAD_transfer(j, 1)
#             Next j
#             Worksheets("Irrigation_summary_" + Irrigation_number).Range("A1").Formula = "#_of_rows"
#             Worksheets("Irrigation_summary_" + Irrigation_number).Range("B1").Formula = Number_of_rows
#         End If
#     Next i
#
#     For i = 1 To Num_treatments
#         If Perform_analysis(i) = 1 Then
#             If i < 10 Then
#                 Irrigation_number = "0" + CStr(i)
#             Else
#                 Irrigation_number = CStr(i)
#             End If
#             Num_days = Worksheets("Irrigation_summary_" + Irrigation_number).Range("B1").Formula
#             Worksheets("Irrigation_summary_" + Irrigation_number).Activate
#
#             For j = 1 To Num_days
#                 Total_depletion = 0
#                 For k = 1 To Num_Cells
#                     Depletions(k) = Worksheets("Irrigation_summary_" + Irrigation_number).Cells(j + 4, k + 1).value
#                     Total_depletion = Depletions(k) + Total_depletion
#                 Next k
#                 Average_depletion(j) = Total_depletion / Num_Cells
#                 Worksheets("Irrigation_summary_" + Irrigation_number).Cells(j + 4, 18).value = Average_depletion(j)
#                 Tot_P = 0
#                 For k = 1 To Num_Cells
#                     Sort_depletions(k) = 100
#
#                     For kk = 1 To Num_Cells
#                         Already_taken = 0
#                         If Depletions(kk) < Sort_depletions(k) Then
#                             For kkk = 2 To k
#                                 If Low_sequence(kkk - 1) = kk Then
#                                     Already_taken = 1
#                                 End If
#                             Next kkk
#
#                             If Already_taken = 0 Then
#                                 Sort_depletions(k) = Depletions(kk)
#                                 Low_sequence(k) = kk
#                             End If
#                         End If
#                     Next kk
#                     Tot_P = PMAD(i, k, j) + Tot_P
#                     Worksheets("Irrigation_summary_" + Irrigation_number).Cells(j + 4, 27).value = Tot_P / Num_Cells * 100
#                 Next k
#                 Ave_P = Tot_P / Num_Cells * 100  'This will be the final P (last j) so it is used to estimate days below
#
#                 Num_in_quartile = 0
#                 Count = 0
#                 For k = 1 To Num_Cells / 4
#                     Total_in_quartile(k) = 0
#                     Total_mm_in_quartile(k) = 0
#                     Quartile_count = 0
#                     For kk = 1 To Num_Cells / 4
#                         Count = Count + 1
#                         Quartile_count = Quartile_count + 1
#                         Total_in_quartile(k) = Total_in_quartile(k) + Depletions(Low_sequence(Count))
#                         Total_mm_in_quartile(k) = Total_mm_in_quartile(k) + mm_depletion(i, Low_sequence(Count), j)
#                     Next kk
#                     Average_in_quartile(k, j) = Total_in_quartile(k) / Quartile_count
#                     Average_mm_in_quartile(k, j) = Total_mm_in_quartile(k) / Quartile_count
#                     Worksheets("Irrigation_summary_" + Irrigation_number).Cells(j + 4, 18 + k).value = Average_in_quartile(k, j)
#                     Worksheets("Irrigation_summary_" + Irrigation_number).Cells(j + 4, 22 + k).value = Average_mm_in_quartile(k, j)
#                 Next k
#                 If j > 2 Then
#                     Ave_irr_required(j) = (Average_mm_in_quartile(1, j) + Average_mm_in_quartile(2, j) + _
#                         Average_mm_in_quartile(3, j) + Average_mm_in_quartile(4, j)) / 4 / Efficiency
#                     Worksheets("Irrigation_summary_" + Irrigation_number).Cells(j + 4, 31).value = Ave_irr_required(j)
#                 End If
#             Next j
#             Average_irr_today_based_on_low_quartile = Average_mm_in_quartile(1, Num_days)
#             Average_irr_today_based_on_high_quartile = Average_mm_in_quartile(4, Num_days)
#             Average_irr_today_based_on_average = (Average_mm_in_quartile(1, Num_days) + Average_mm_in_quartile(2, Num_days) + _
#             Average_mm_in_quartile(3, Num_days) + Average_mm_in_quartile(4, Num_days)) / 4
#
#             start_day = days(2, 1)
#
#             For k = 1 To 10
#                 If Num_days > 4 Then
#                     rowi = Num_days + 1
#                     rowf = Num_days + 4
#
#                     Select Case k
#                     Case Is = 1
#                         Column = "R"
#                         Column_number = 18
#                     Case Is = 2
#                         Column = "S"
#                         Column_number = 19
#                     Case Is = 3
#                         Column = "T"
#                         Column_number = 20
#                     Case Is = 4
#                         Column = "U"
#                         Column_number = 21
#                     Case Is = 5
#                         Column = "V"
#                         Column_number = 22
#                     Case Is = 6
#                         Column = "W"
#                         Column_number = 23
#                     Case Is = 7
#                         Column = "X"
#                         Column_number = 24
#                     Case Is = 8
#                         Column = "Y"
#                         Column_number = 25
#                     Case Is = 9
#                         Column = "Z"
#                         Column_number = 26
#                     Case Is = 10
#                         Column = "AA"
#                         Column_number = 27
#                     End Select
#                     For kk = 1 To 4
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Cells(rowf + kk, 17).value = kk
#                     Next kk
#                     x_address = Column + CStr(rowi) + ":" + Column + CStr(rowf)
#                     rx = Range(x_address)
#                     y_address = "Q" + CStr(rowi + 4) + ":Q" + CStr(rowf + 4)
#                     ry = Range(y_address)
#                     slope(k) = WorksheetFunction.slope(rx, ry)
#                     yint(k) = WorksheetFunction.Intercept(rx, ry)
#                     For kk = 1 To 4
#                         Projected_depletion = (kk + 4) * slope(k) + yint(k)
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Cells(rowf + kk, Column_number).value = Projected_depletion
#                     Next kk
#                     If k = 1 Then
#                         Present_percent_dep = Worksheets("Irrigation_summary_" + Irrigation_number).Cells(rowf, Column_number).value
#                         Days_to_irrigate(k) = (Ave_P - Present_percent_dep) / slope(k)    'P  = present + slope * days    --> Days = (P - present) / slope
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AH4").value = CInt(Days_to_irrigate(k))
#                         Projected_DOY = Current_day_of_year + CInt(Days_to_irrigate(k))
#                         If Projected_DOY > 365 Then
#                             Projected_DOY = Projected_DOY - 365
#                         End If
#                         Place_date = Worksheets("DOY dates").Cells(Projected_DOY + 2, 2).value
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AI4").value = Place_date
#                     ElseIf k = 2 Then
#                         Present_percent_dep = Worksheets("Irrigation_summary_" + Irrigation_number).Cells(rowf, Column_number).value
#                         Days_to_irrigate(k) = (Ave_P - Present_percent_dep) / slope(k)    'P  = present + slope * days    --> Days = (P - present) / slope
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AE4").value = CInt(Days_to_irrigate(k))
#                         Projected_DOY = Current_day_of_year + CInt(Days_to_irrigate(k))
#                         If Projected_DOY > 365 Then
#                             Projected_DOY = Projected_DOY - 365
#                         End If
#                         Place_date = Worksheets("DOY dates").Cells(Projected_DOY + 2, 2).value
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AF4").value = Place_date
#                     ElseIf k = 5 Then
#                         Present_percent_dep = Worksheets("Irrigation_summary_" + Irrigation_number).Cells(rowf, Column_number).value
#                         Days_to_irrigate(k) = (Ave_P - Present_percent_dep) / slope(k)    'P  = present + slope * days    --> Days = (P - present) / slope
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AB4").value = CInt(Days_to_irrigate(k))
#                         Projected_DOY = Current_day_of_year + CInt(Days_to_irrigate(k))
#                         If Projected_DOY > 365 Then
#                             Projected_DOY = Projected_DOY - 365
#                         End If
#                         Place_date = Worksheets("DOY dates").Cells(Projected_DOY + 2, 2).value
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AC4").value = Place_date
#                     ElseIf k = 6 Then
#                         Present_depth_low = Worksheets("Irrigation_summary_" + Irrigation_number).Cells(rowf, Column_number).value
#                         Depth_to_irrigate = Present_depth_low + slope(k) * Days_to_irrigate(2)
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AG4").value = Depth_to_irrigate
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AG5").value = Average_irr_today_based_on_low_quartile
#                     ElseIf k = 9 Then
#                         Present_depth_high = Worksheets("Irrigation_summary_" + Irrigation_number).Cells(rowf, Column_number).value
#                         Depth_to_irrigate = Present_depth_high + slope(k) * Days_to_irrigate(5)
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AD4").value = Depth_to_irrigate
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AD5").value = Average_irr_today_based_on_high_quartile
#                     ElseIf k = 10 Then
#                         Present_depth = (Present_depth_low + Present_depth_high) / 2
#                         Depth_to_irrigate = Present_depth + (slope(6) + slope(9)) / 2 * Days_to_irrigate(1)
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AJ4").value = Depth_to_irrigate
#                         Worksheets("Irrigation_summary_" + Irrigation_number).Range("AJ5").value = Average_irr_today_based_on_average
#                     End If
#
#                 End If
#
#             Next k
#         End If
#     Next i
#
#     self.Plot_analysis As Boolean
#     Plot_analysis = 0
#     For i = 1 To Num_plots
#         Perform_plot_analysis(i) = Worksheets("Spatial_data").Cells(i + 26, 13).value
#         If Perform_plot_analysis(i) = 1 Then
#             Plot_analysis = 1
#         End If
#     Next i
#
#     If Plot_analysis = 1 Then
#         Worksheets("Plot irrigation summaries").Range("A5:A1500").ClearContents
#         Worksheets("Plot irrigation summaries").Range("C5:E1500").ClearContents
#     End If
#     For i = 1 To Num_plots
#         If Perform_plot_analysis(i) = 1 Then
#             Plot_number = Worksheets("Spatial_data").Cells(26 + i, 12)
#             If Planting_date_in_spatial_data Then
#                 Planting_day = Worksheets("Spatial_data").Cells(26 + i, 16).value
#             Else
#                 Planting_day = Crop_data_planting_days(1, Plot_number + 1)
#             End If
#             Current_day_of_year = Worksheets("Main").Range("E3").value
#             Number_of_rows = Current_day_of_year - Planting_day + 2
#             Cell_range_percent_depletion = "IR3:IR" + CStr(Number_of_rows + 3)
#             Cell_range_percent_depletion_out = "D5:D" + CStr(Number_of_rows + 5)
#             Cell_range_mm_depletion = "IQ3:IQ" + CStr(Number_of_rows + 3)
#             Cell_range_mm_depletion_out = "C5:C" + CStr(Number_of_rows + 5)
#             Cell_range_days = "W3:W" + CStr(Number_of_rows + 3)
#             Cell_range_days_out = "A5:A" + CStr(Number_of_rows + 4)
#             Cell_range_pmad = "IT3:IT" + CStr(Number_of_rows + 3)
#             Cell_range_pmad_out = "E5:E" + CStr(Number_of_rows + 4)
#             Irrigation_name = Worksheets("Spatial_data").Cells(26 + Plot_number, 17).Formula
#             Irrigation_depths = Worksheets("Irr_" + Irrigation_number).Range("C1:C1200")
#             Irrigation_number = Right(Irrigation_name, 2)
#             Section_number = Worksheets("Spatial_data").Cells(26 + Plot_number, 15).value
#             Worksheet_name = "C_" + CStr(Plot_number)
#             Worksheets("Plot irrigation summaries").Range("A3").Formula = Worksheet_name
#             days = Worksheets(Worksheet_name).Range(Cell_range_days)
#             days_right = Worksheets(Worksheet_name).Range(Cell_range_days_right)
#             Percent_depletion_transfer = Worksheets(Worksheet_name).Range(Cell_range_percent_depletion)
#             PMAD_transfer = Worksheets(Worksheet_name).Range(Cell_range_pmad)
#             mm_depletion_transfer = Worksheets(Worksheet_name).Range(Cell_range_mm_depletion)
#             Worksheets("Plot irrigation summaries").Range(Cell_range_days_out) = days
#             Worksheets("Plot irrigation summaries").Range(Cell_range_mm_depletion_out) = mm_depletion_transfer
#             Worksheets("Plot irrigation summaries").Range(Cell_range_percent_depletion_out) = Percent_depletion_transfer
#             Worksheets("Plot irrigation summaries").Range(Cell_range_pmad_out) = PMAD_transfer
#             Efficiency_one_plot = Worksheets("Plot irrigation summaries").Range("D1")
#
#             Worksheets("Plot irrigation summaries").Range("A1").Formula = "#_of_rows"
#             Worksheets("Plot irrigation summaries").Range("B1").Formula = Number_of_rows
#             Num_days = Worksheets("Plot irrigation summaries").Range("B1").value
#             Column = "A"
#             Column_number = 4
#
#             If Num_days > 4 Then
#                 rowi = Num_days + 1
#                 rowf = Num_days + 4
#                 Final_mm = Worksheets("Plot irrigation summaries").Cells(rowf, 3).value
#                 Final_percent = Worksheets("Plot irrigation summaries").Cells(rowf, 4).value
#                 Ratio = Final_mm / Final_percent
#                 For kk = 1 To 4
#                     Worksheets("Plot irrigation summaries").Cells(rowf + kk, 1).value = kk
#                     Worksheets("Plot irrigation summaries").Cells(rowf + kk, 2).value = Worksheets("Plot irrigation summaries").Cells(rowf + kk - 1, 2).value + 1
#                 Next kk
#                 x_address = Column + CStr(rowi + 4) + ":" + Column + CStr(rowf + 4)
#                 rx = Range(x_address)
#                 y_address = "D" + CStr(rowi) + ":D" + CStr(rowf)
#                 ry = Range(y_address)
#                 slope(1) = WorksheetFunction.slope(ry, rx)
#                 yint(1) = WorksheetFunction.Intercept(ry, rx)
#                 For kk = 1 To 50
#                     Projected_depletion = (kk + 4) * slope(1) + yint(1)
#                     Worksheets("Plot irrigation summaries").Cells(rowf + kk, Column_number).value = Projected_depletion
#                     Projected_mm_depletion = Projected_depletion * Ratio
#                     Worksheets("Plot irrigation summaries").Cells(rowf + kk, Column_number - 1).value = Projected_mm_depletion
#                 Next kk
#                 Days_to_irrigate_one_plot = -(Worksheets("Plot irrigation summaries").Cells(rowf, 4).value - _
#                     Worksheets("Plot irrigation summaries").Cells(rowf, 5).value * 100) / slope(1) 'P  = present + slope * days    --> Days = (P - present) / slope
#                 Worksheets("Plot irrigation summaries").Range("B2").value = CInt(Days_to_irrigate_one_plot)
#                 Place_date = Worksheets("Plot irrigation summaries").Cells(rowf + CInt(Days_to_irrigate_one_plot), 2).value
#                 Worksheets("Plot irrigation summaries").Range("D2").value = Place_date
#                 Projected_depletion = Worksheets("Plot irrigation summaries").Cells(rowf + CInt(Days_to_irrigate_one_plot), 3).value
#                 Worksheets("Plot irrigation summaries").Range("F2").value = Projected_depletion / Efficiency_one_plot
#                 Depth_today = Worksheets("Plot irrigation summaries").Cells(rowf, 3) / Efficiency_one_plot
#                 Worksheets("Plot irrigation summaries").Range("F1").value = Depth_today
#                 Clear_zone = "A" + CStr(rowf + 1) + ":A" + CStr(rowf + 4)
#                 Worksheets("Plot irrigation summaries").Range(Clear_zone).ClearContents
#             End If
#         End If
#     Next i
# End Sub