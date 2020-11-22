import requests
import urllib.parse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from math import *

frame_dict = {1: [0.028,0.8], 2: [0.016, 0.9], 3: [0.03, 0.75], 4: [0.02, 0.75]}

def main():
    print('Welcome to EQ load design program for symmetrical buildings')
    usgs_result_dict = read_earthquake_usgs_gov()
    print(usgs_result_dict)
    soil_t = read_site_class()
    s1_coef = initial_s_coefficients('max_direction_s1', 'mapped_s1', 'deterministic_floor_s1', 'cr1')
    ss_coef = initial_s_coefficients('max_direction_ss', 'mapped_ss', 'deterministic_floor_ss', 'crs')
    sds = (2/3)*ss_coef*soil_adjustment_coefficient_fa(ss_coef, soil_t)
    sd1 = (2/3)*s1_coef*soil_adjustment_coefficient_f1(s1_coef, soil_t)
    list_of_period_values = period_values(sd1,sds,usgs_result_dict)
    list_of_acceleration_values = acceleration_values(sd1,sds,list_of_period_values)
    seismic_response_map = plot_response_spectrum(list_of_period_values,list_of_acceleration_values )
    
    building_length = float(input("Enter length of the building in feet "))
    # add a note, what is defined as length
    building_width = float(input("Enter width of the building in feet "))
    bay_width = float(input("Enter width of the bay in feet "))
    number_of_storeys = int(input("Enter number of storeys "))
    height_of_storeys = float(input("Enter height of the storey "))
    number_of_braces_eq_dir = int(input("Enter number of braces in EQ direction "))
    number_of_braces_noneq_dir = int(input("Enter number of braces in nonEQ direction "))
    slope_roof = int(input("Enter the value for slope of the roof "))
    dead_loads = int(input("Enter the DL value (per floor) in psf "))
    dead_loads_roof = int(input("Enter the roof DL value in psf "))
    live_loads = int(input("Enter the LL value (per floor) in psf "))
    live_loads_roof = int(input("Enter the LL_r(roof) value (per floor) in psf "))
    moment_frame = int(input("Choose the type of moment-resisting frame (1 - steel, 2 - concrete, 3 - eccentrically braced, 4 - other structural system"))
    r_coeff = int(input("Enter Response modification coefficient from the Table 12.2.1 of ASCE "))
    importance = int(input("Enter Response Importance factor from the Table 1.5.2 of ASCE "))
    distance_braces_eq = int(input("Enter distance between braces in EQ direction "))
    distance_non_eq = int(input("Enter distance between braces in non-EQ direction "))
    torsion_coeff = 0.05*building_length
    total_k = number_of_braces_eq_dir*((distance_braces_eq/2)**2) + number_of_braces_noneq_dir*((distance_non_eq/2)**2)
    area = building_length * building_width
    weight = total_weight(area,dead_loads,dead_loads_roof, number_of_storeys)
    period_of_building = frame_dict[moment_frame][0]*((number_of_storeys*height_of_storeys)**frame_dict[moment_frame][1])
    print("Period of the building is",period_of_building, "sec")
    cs_coeff = calculate_cs_coefficient(period_of_building)
    base_shear = cs_coeff*weight
    print("Base shear of the building is",base_shear, "kip")
    weigh_h_dict = calculate_base_shear_distribution(period_of_building)
    lateral_distribution_coeff = calculate_lateral_distribution_coefficient(period_of_building,base_shear)
    print("")
    print("This part of the program computes lateral force distribution")
    print("Lateral distribution coeefficients (Cvx) are:",lateral_distribution_coeff)
    floor_lateral_force = calculate_floor_lateral_force(lateral_distribution_coeff, base_shear)
    print("Lateral force distribution per floor:",floor_lateral_force)
    brace_lateral_force = calculate_lateral_force_brace(floor_lateral_force)
    print("Lateral force distribution per brace:",brace_lateral_force)
    print("Accidental torsion for the building: T = " + torsion_coeff + "F")
    print("Total stiffness of the structure: k_t = " + total_k + "K")
    f_eq_dir = ((distance_braces_eq/2)*torsion_coeff)/total_k
    f_noneq_dir = ((distance_non_eq/2)*torsion_coeff)/total_k
    print("Lateral force due to torsion in the EQ direction (per floor): Fx = " + f_eq_dir + "Fi")
    print("Lateral force due to torsion in the EQ direction (per floor): Fx = " + f_noneq_dir + "Fi")
    torsion_brace_lateral_force = calulate_torsion_lateral_force_distribution(brace_lateral_force)
    print("Lateral force distribution per brace with account for accidental torsion:", torsion_brace_lateral_force)
    # Vertical loads
    total_dead = ((dead_loads*number_of_storeys + dead_loads_roof)*(bay_width**2))/1000
    print("Total force from DL is",total_dead, "kip")
    ll_force = calculate_ll()
    ll_force_roof = calculate_ll_roof()
    ll_force_total = ll_force + ll_force_roof
    print("Force from LL is",ll_force, "kip")
    print("Force from roof LL is",ll_force_roof, "kip")
    print("Total force from LL is",ll_force_total, "kip")

def read_earthquake_usgs_gov():
    longitude, latitude = read_coordinates()
    reference_document = read_reference_document()
    site_class = read_site_class()
    risk_category = read_risk_category()

    r = requests.get('https://earthquake.usgs.gov/designmaps/beta/us/service/' 
        + reference_document + '/' + site_class + '/' + risk_category + '/' + longitude + '/' + latitude + '/t')
    return r.json()

def read_coordinates():
    print('(1/4) Please provide either')
    print('  1: a street address')
    print('  2: longitude and latitude')
    
    geo_type = read_geo_type()
    if geo_type == '1':
        address = input('Address: ')
        return find_coordinates(address)
    elif geo_type == '2':
        return read_longitude(), read_latitude()

def read_geo_type():
    while True:
        geo_type = input('?: ')
        if is_geo_type_valid(geo_type):
            return geo_type    
        else:
            print('ERROR: invalid type of geo location')

def is_geo_type_valid(geo_type):
    # TODO: implement!
    return True;

def find_coordinates(address):
    r = requests.get('https://geocode.arcgis.com/arcgis/rest/services/World/GeocodeServer/find?f=json&text=' + urllib.parse.quote_plus(address))
    response = r.json()
    return str(response['locations'][0]['feature']['geometry']['x']), str(response['locations'][0]['feature']['geometry']['y'])

def read_longitude():
    while True:
        longitude = input('Longitude: ')
        if is_longitude_valid(longitude):
            return longitude
        else:
            print('ERROR: Invalid longitude value')    

def is_longitude_valid(longitude):
    # TODO: implement!
    return True;

def read_latitude():
    while True:
        latitude = input('Latitude: ')
        if is_latitude_valid(latitude):
            break
        else:
            print('ERROR: Invalid latitude value')
    return latitude    

def is_latitude_valid(latitude):
    # TODO: implement!
    return True;

def read_reference_document():
    print('(2/4) We will be using "2015 NEHRP Provisions" as a reference document')
    return '1'

def read_site_class():
    print('(3/4) There are 7 site classes to choose from:')
    print('  1: Hard Rock')
    print('  2: Rock (measured)')
    print('  3: Rock (unmeasured)')
    print('  4: Very Dense Soil and Soft Rock')
    print('  5: Stiff Soil (determined)')
    print('  6: Stiff Soil (default)')
    print('  7: Soft Clay Soil')
    while True:
        site_class = input('?: ')
        if is_site_class_valid(site_class):
            break
        else:
            print('ERROR: invalid site class')
    return site_class

def is_site_class_valid(site_class):
    # TODO: implement!
    return True;

def read_risk_category():
    print('(4/4) There are 4 risk categories to choose from:')
    print('  1: I')
    print('  2: II')
    print('  3: III')
    print('  4: IV (Essential Facilities)')
    while True:
        risk_category = input('?: ')
        if is_risk_category_valid(risk_category):
            break
        else:
            print('ERROR: invalid risk category')
    return risk_category

def is_risk_category_valid(risk_category):
    # TODO: implement!
    return True;

def s_avg_value(name, dictionary):
    lst = []
    for el in dictionary:
        lst.append(el[name])
    return np.average(lst)

def initial_s_coefficients(max_direction, mapped_val, deterministic_floor_val, cr_val):
    s_val = 0
    dir_coef_s = usgs_result_dict['output']['metadata'][max_direction]
    output_dict = usgs_result_dict['output']['data']
    s_var = s_avg_value(mapped_val, output_dict)
    s_uh = dir_coef_s*s_var
    cr = s_avg_value(cr_val, output_dict)
    floor_s = usgs_result_dict['output']['metadata'][deterministic_floor_val]
    s_val = min(s_uh*cr, floor_s)
    return s_val  

def linear_interpolation(x1,x2,y1,y2,val):
    xp = [x1,x2]
    fp = [y1,y2]
    return np.interp(val, xp, fp)

def soil_adjustment_coefficient_fa(ss, soil_type):
    # pass list containing Ss and S1, and soil_type
    f_a = 0
    if soil_type == 1:
        f_a = 0.8
    elif soil_type == 2:
        f_a = 0.9
    elif soil_type == 3:
        f_a = 1.0
    elif soil_type == 4:
        if ss <=0.5:
            f_a = 1.3
        elif ss >= 0.75:
            f_a = 1.2
        else:
            f_a = linear_interpolation(0.5,0.75,1.3,1.2,ss)
    elif soil_type == 5:
        if ss <= 0.25:
            f_a = 1.6
        elif ss == 0.5:
            f_a = 1.4
        elif 0.25 < ss < 0.5:
            f_a = linear_interpolation(0.25,0.5,1.6,1.4, ss)
        elif ss == 0.75:
            f_a = 1.2
        elif 0.5 < ss < 1:
            f_a = linear_interpolation(0.5,0.75,1.4,1.2, ss)
        elif ss == 1.0:
            f_a = 1.1
        elif 0.75 < ss < 1:
            f_a = linear_interpolation(0.75,1.0,1.2,1.1, ss)
        elif ss <= 1.25:
            f_a = 1.0
        else:
            f_a = linear_interpolation(1.0,1.25,1.1,1.0, ss)
    elif soil_type == 6:
        if ss <= 0.25:
            f_a = 1.6
        elif ss == 0.5:
            f_a = 1.4
        elif 0.25 < ss < 0.5:
            f_a = linear_interpolation(0.25,0.5,1.6,1.4, ss)
        elif ss <= 0.75:
            f_a = 1.2
        else:
            f_a = linear_interpolation(0.5,0.75,1.4,1.2, ss)
    else:
        if ss <= 0.25:
            f_a = 2.4
        elif ss == 0.5:
            f_a = 1.7
        elif 0.25 < ss < 0.5:
            f_a = linear_interpolation(0.25,0.5,2.4,1.7, ss)
        elif ss == 0.75:
            f_a = 1.3
        elif 0.5 < ss < 0.75:
            f_a = linear_interpolation(0.5,0.75,1.7,1.3, ss)
        elif ss <= 1.0:
            f_a = 1.2
        else:
            f_a = linear_interpolation(0.75,1.0,1.3,1.2, ss)
    return f_a

def soil_adjustment_coefficient_f1(s1, soil_type):
    f_1 = 0
    if soil_type == 1:
        f_1 = 0.8
    elif soil_type == 2:
        f_1 = 0.9
    elif soil_type == 3:
        f_1 = 1.0
    elif soil_type == 4:
        if s1 <= 0.5:
            f_1 = 1.5
        elif s1 >= 0.6:
            f_1 = 1.4
        else:
            f_1 = linear_interpolation(0.5,0.6,1.5,1.4,s1)
    elif soil_type == 5 or soil_type == 6:
        if s1 <= 0.1:
            f_1 = 2.4
        elif s1 == 0.2:
            f_1 = 2.2
        elif s1 > 0.1 and s1 < 0.2:
            f_1 = linear_interpolation(0.1,0.2,2.4,2.2,s1)
        elif s1 == 0.3:
            f_1 = 2.0
        elif s1 == 0.4:
            f_1 = 1.9
        elif s1 > 0.3 and s1 < 0.4:
            f_1 = linear_interpolation(0.3,0.4,2.0,1.9,s1)
        elif s1 == 0.5:
            f_1 = 1.8
        elif s1 >= 0.6:
            f_1 = 1.7
        else:
            f_1 = linear_interpolation(0.5,0.6,1.8,1.7,s1)
    else:
        if s1 <= 0.1:
            f_1 = 4.2
        elif s1 == 0.2:
            f_1 = 3.3
        elif s1 > 0.1 and s1 < 0.2:
            f_1 = linear_interpolation(0.1,0.2,4.2,3.3,s1)
        elif s1 == 0.3:
            f_1 = 2.8
        elif s1 == 0.4:
            f_1 = 2.4
        elif s1 > 0.3 and s1 < 0.4:
            f_1 = linear_interpolation(0.3,0.4,2.8,2.4,s1)
        elif s1 == 0.5:
            f_1 = 2.2
        elif s1 >= 0.6:
            f_1 = 2.0
        else:
            f_1 = linear_interpolation(0.5,0.6,2.2,2.0,s1)
    return f_1
        
def period_values(sd1,sds,result_dict):
    t_zero = 0.2*(sd1/sds)
    t_short =(sd1/sds)
    t_long = usgs_result_dict['output']['tl']
    period = [t_zero,t_short,t_long]
    return period

def acceleration_values(sd1,sds,period_list):
    s_initial = sds*0.4
    s_t0 = sds
    s_ts = sds
    s_tl = sd1/period_list[3]
    acceleration = [s_initial,s_t0,s_ts,s_tl]
    return acceleration

def plot_response_spectrum(pl, accel):
    plt.clf()
    xs = pl
    ys = accel
    plt.plot(xs, ys)
    plt.xscale('symlog')
    plt.title("Acceleration response spectrum")
    plt.xlim(-0.05,6)
    plt.xlabel('Period, seconds')
    plt.ylabel('Acceleration (g)')
    plt.show()

def total_weight(area,dl,dlr,storeys):
    storey_weigth = float(area*dl) / 1000
    roof_weight = float(area*dlr) / 1000
    return storey_weigth*storeys + roof_weight

def calculate_cs_coefficient(period):
    c_s = sds/(r_coef/importance)
    if c_s < 0.1:
        c_s = 0.1
    elif c_s > 0.6:
        c_s = 0.6
    elif period < period_list[3]:
        if c_s > sds/(period*(r_coef/importance)):
            c_s = sds/(period*(r_coef/importance))
    elif period > period_list[3]:
        if c_s > (sds*period)/((period**2)*(r_coef/importance)):
            c_s = (sds*period)/((period**2)*(r_coef/importance))
    return c_s

def calculate_base_shear_distribution(period):
    if period <=0.5:
        k = 1
    elif period >= 2.5:
        k = 2
    else:
        k = linear_interpolation(0.5,2.5,1,2,period)
    weight_height_numerator = {}
    total_val = 0
    floor_val = 0
    roof_val = 0
    for floor in range(number_of_storeys):
        floor_val += (dead_loads*area*(height_of_storeys**k))/1000
        weight_height_numerator[floor] = floor_val
    roof_val = (dead_loads_roof*area*(height_of_storeys**k)*(number_of_storeys +1))/1000
    weight_height_numerator.update({'roof' : roof_val})
    return weight_height_numerator 

def calculate_lateral_distribution_coefficient(period,shear):
    #weigh_h_dict = calculate_base_shear_distribution()
    #print(sum(weigh_h_dict.values()))
    values = weigh_h_dict.values()
    keys = weigh_h_dict.keys()
    lateral_coef_dict = { k:v/(sum(weigh_h_dict.values())) for (k,v) in zip(keys, values)}
    return lateral_coef_dict 

def calculate_floor_lateral_force(coefficient, shear):
    lateral_force_dict = {key:value*shear for (key,value) in coefficient.items()}
    return lateral_force_dict

def calculate_lateral_force_brace(lateral_f):
    lateral_force_brace_dict = {key:value/number_of_braces_eq_dir for (key,value) in lateral_f.items()}
    return lateral_force_brace_dict

def calulate_torsion_lateral_force_distribution(b_lat_force):
    tor_lateral_force_brace_dict = {key:(value + f_noneq_dir/2) for (key,value) in b_lat_force.items()}
    return tor_lateral_force_brace_dict

def calculate_ll():
    if (bay_width**2)*number_of_storeys < 400:
        ll_force = live_loads*(bay_width**2)*number_of_storeys
    else:
        red_factor = (0.25 + (15 / math.sqrt(4*(bay_width**2)*number_of_storeys)))
        ll_force = ((bay_width**2)*number_of_storeys*red_factor*live_loads)/1000
    return ll_force

def calculate_ll_roof():
    r1 = 0
    r2 = 0
    if (bay_width**2) <= 200:
        r1 = 1
    elif (bay_width**2) >= 600:
        r1 = 0.6
    else:
        r1 = 1.2 - 0.001*(bay_width**2)
    if slope_roof <=4:
        r2 = 1
    elif slope_roof >= 12:
        r2 = 0.6
    else:
        r2 = 1.2 - 0.05*slope_roof
    ll_force_r = (r1*r2*live_loads_roof*(bay_width**2))/1000
    return ll_force_r


def base_shear():
    shear_value = total_load * seismic response coefficient
    return shear_value


if __name__ == "__main__":
    main()