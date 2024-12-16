


#11/14 updates mean after each cycle and not just comapred to pert

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pickle
#Load the dictionary from the pickle file
with open('data_stream_AB03', 'rb') as f:
    loaded_data = pickle.load(f)


false_det = 0
no_det = 0

pert = 2000
compared_cyc = []
# Access the variables
time = loaded_data['time']
COM_ALL = loaded_data['COM_ALL']
RHL_ALL = loaded_data['RHL_ALL']
LHL_ALL = loaded_data['LHL_ALL']
VEL_COM_ALL = loaded_data['VEL_COM_ALL']
VEL_RHL_ALL = loaded_data['VEL_RHL_ALL']
VEL_LHL_ALL = loaded_data['VEL_LHL_ALL']
indxr = loaded_data['indxr']
psum_cont = []

last_force = 0
heel_strike = False
frame = 0
heel_strike_count = 0
cont_c_rel = []
cont_r_rel = []
cont_l_rel = []
#raw cont
cont_c_raw = []
cont_r_raw = []
cont_l_raw = []
#ony current cycle values
cyc_c = []
cyc_r = []
cyc_l = []

#raw current cycle
cyc_c_raw = []
cyc_r_raw = []
cyc_l_raw = []
cyc_len_list = []
#split cycles, cont
c_s = []
r_s = []
l_s = []
#velocity over each cycle 
c_s_v = []
r_s_v =[]
l_s_v = []
#only this cycle
cyc_cv = []
cyc_rv = []
cyc_lv = []
X = []


getted_frames_cycle = []
pert_trigger = False
tracking_active = False
weight_sum = False
early_flag = 0
def A_value(cpmean,cpstd,cpoint):
    A = 0
    lambdaa= 0
    if (cpmean-2*cpstd)>cpoint:
        #cpstd = max(cpstd,0.001)
        A = abs(cpoint-(cpmean-2*cpstd)) 
        lambdaa = (1/16)*(1/(2*abs(cpstd)+A))
    if (cpmean+2*cpstd)<cpoint:
        #cpstd = max(cpstd,0.001)
        A = abs(cpoint-(cpmean+2*cpstd))
        lambdaa = (1/16)*(1/(2*abs(cpstd)+A))
    Pm = np.multiply(A,lambdaa)

    return Pm

def A_valuevel(cpmean,cpstd,cpoint):
    A = 0
    lambdaa= 0
    if (cpmean-2*cpstd)>cpoint:
        A = abs(cpoint-(cpmean-2*cpstd)) 
        cpstd = cpstd/cpmean
        lambdaa = (1/16)*(1/(2*abs(cpstd)+A))
    
    if (cpmean+2*cpstd)<cpoint:
        A = abs(cpoint-(cpmean+2*cpstd))
        cpstd = cpstd/cpmean
        lambdaa = (1/16)*(1/(2*abs(cpstd)+A))
    Pm = np.multiply(A,lambdaa)

    return Pm

def resample(positions, num_points=100):
    """
    Resample foot position data to be evenly spaced between 0 and 100% of the gait cycle.

    Parameters:
    foot_positions (list or np.ndarray): A 2D list or array of shape (N, 3), where N is the number of original timestamps,
                                         and each row represents the (x, y, z) position of the foot.
    num_points (int): The number of points to interpolate to, representing the 0–100% gait cycle.

    Returns:
    np.ndarray: A 2D array of shape (num_points, 3) containing the resampled foot positions.
    """
    positions = np.array(positions)  # Ensure input is an array
    num_original_points = positions.shape[0]
   
    # Original timestamps (normalized from 0 to 100% based on the number of input points)
    original_timestamps = np.linspace(0, 100, num_original_points)
   
    # Target timestamps for interpolation (0% to 100% in the specified number of points)
    target_timestamps = np.linspace(0, 100, num_points)
   
    # Interpolate each axis (x, y, z) independently
    x_interp = interp1d(original_timestamps, positions[:, 0])
    y_interp = interp1d(original_timestamps, positions[:, 1])
    z_interp = interp1d(original_timestamps, positions[:, 2])
   
    # Generate resampled points for each axis at the target timestamps
    resampled_x = x_interp(target_timestamps)
    resampled_y = y_interp(target_timestamps)
    resampled_z = z_interp(target_timestamps)
   
    # Combine into a 2D array for the full resampled cycle
    resampled_positions = np.vstack((resampled_x, resampled_y, resampled_z)).T
   
    return resampled_positions

trials = len(COM_ALL)
avgd = []
good_det = 0
false_det = 0
pertflag = 0
for trial in range(trials):
    #plt.plot(list(range(100)),np.swapaxes(np.array(r_s_v)[:,:,1],0,1),color="blue");plt.plot(mean_rhlv[:,1]+2*std_rhlv[:,1],color = 'black');plt.plot(mean_rhlv[:,1]-2*std_rhlv[:,1],color = 'black');plt.plot(mean_rhlv[:,1],color = 'black');plt.plot(list(range(100)),np.array(up_cyc_rv)[:,1],color="red");plt.plot(list(range(100)),np.array(r_s_v)[11,:,1],color="green")
    heel_strike_count = 0
    c_s = []
    r_s = []
    l_s = []
    c_s_v = []
    r_s_v = []
    l_s_v = []
    pind = 0
    delay = 0
    early_flag = 0
    # cyc_c = []
    # cyc_l = []
    # cyc_r = []

    # cyc_cv = []
    # cyc_lv = []
    # cyc_rv = []
    timet0 = time[trial]
    timet = timet0 - timet0[0]
    indxrt = indxr[trial]
    comt = COM_ALL[trial]
    rhlt = RHL_ALL[trial]
    lhlt = LHL_ALL[trial]
    comvelt = VEL_COM_ALL[trial]
    rhlvelt = VEL_RHL_ALL[trial]
    lhlvelt = VEL_LHL_ALL[trial]
    for ind,t in enumerate(timet):
        if pert == 1:
            break
        com_rel = comt[ind,:]
        rhl_rel = rhlt[ind,:]
        lhl_rel = lhlt[ind,:]
        if t > timet[0] and ind<len(comt):
            cur_v_c = comvelt[ind,:]
            cur_v_r = rhlvelt[ind,:]
            cur_v_l = lhlvelt[ind,:]
            if heel_strike_count >= len(indxrt):
                    break
            heel_strike = 1 if t == indxrt[heel_strike_count] else 0
            if heel_strike:
                heel_strike_count +=1
                if heel_strike_count >= len(indxrt):
                    break
                print("strike")
                if tracking_active:
                    cyc_len,_ = np.shape(cyc_c)
                    cyc_len_list.append(cyc_len)
                    #mean and std first 10 cycles
                    if heel_strike_count == 4:
                        #COM
                        up_cyc_c =resample(np.array(cyc_c))
                        c_s.append(up_cyc_c)
                        #plt.plot(list(range(100)),np.swapaxes(np.array(c_s)[:,:,0],0,1))
                        mean_com = up_cyc_c
                        #RIGHT
                        up_cyc_r =resample(np.array(cyc_r))
                        r_s.append(up_cyc_r)
                        mean_rhl = up_cyc_r
                        #LEFT
                        up_cyc_l =resample(np.array(cyc_l))
                        l_s.append(up_cyc_l)
                        mean_lhl = up_cyc_l
                        #append to list for making a continous mean 

                        #Velocity 
                        up_cyc_cv =resample(cyc_cv)
                        #add to spi cycles vector for later mean
                        c_s_v.append(up_cyc_cv)
                        mean_comv = up_cyc_cv

                        up_cyc_rv =resample(cyc_rv)
                        r_s_v.append(up_cyc_rv)
                        mean_rhlv = up_cyc_rv 

                        up_cyc_lv =resample(cyc_lv)
                        l_s_v.append(up_cyc_lv)
                        mean_lhlv = up_cyc_lv 
                        
                    if ( heel_strike_count >4  and heel_strike_count<=14):
                        up_cyc_c = resample(np.array(cyc_c))
                        c_s.append(up_cyc_c)
                        mean_com = np.mean((np.array(c_s)),axis=0)
                        std_com = np.std((np.array(c_s)),axis=0)

                        up_cyc_r = resample(np.array(cyc_r))
                        r_s.append(up_cyc_r)
                        mean_rhl = np.mean((np.array(r_s)),axis=0)
                        std_rhl = np.std((np.array(r_s)),axis=0)

                        up_cyc_l = resample(np.array(cyc_l))
                        l_s.append(up_cyc_l)
                        mean_lhl = np.mean((np.array(l_s)),axis=0)
                        std_lhl = np.std((np.array(l_s)),axis=0)
                        
                        #Velocity 
                        up_cyc_cv =resample(cyc_cv)
                        #add to spi cycles vector for later mean
                        c_s_v.append(up_cyc_cv)
                        mean_comv = np.mean((np.array(c_s_v)),axis=0)
                        std_comv = np.std((np.array(c_s_v)),axis=0)

                        up_cyc_rv =resample(cyc_rv)
                        r_s_v.append(up_cyc_rv)
                        mean_rhlv = np.mean((np.array(r_s_v)),axis=0)
                        std_rhlv = np.std((np.array(r_s_v)),axis=0)

                        up_cyc_lv =resample(cyc_lv)
                        l_s_v.append(up_cyc_lv)
                        mean_lhlv = np.mean((np.array(l_s_v)),axis=0)
                        std_lhlv = np.std((np.array(l_s_v)),axis=0)

                    if heel_strike_count>14:
                        # if early_flag ==1:
                        #     break
                        weight_sum = True
                        up_cyc_c = resample(np.array(cyc_c))
                        c_s.append(up_cyc_c)
                        
                        up_cyc_r = resample(np.array(cyc_r))
                        r_s.append(up_cyc_r)
                        
                        up_cyc_l = resample(np.array(cyc_l))
                        l_s.append(up_cyc_l)
                        #Velocity 
                        up_cyc_cv =resample(cyc_cv)
                        #add to spi cycles vector for later mean
                        c_s_v.append(up_cyc_cv)
                        up_cyc_rv =resample(cyc_rv)
                        r_s_v.append(up_cyc_rv)
                        
                        up_cyc_lv =resample(cyc_lv)
                        l_s_v.append(up_cyc_lv)
                        psum_cont = []
                        for compare_index in range (100):
                            # if early_flag ==1:
                            #     break
                            PCX = A_valuevel(mean_com[compare_index,0],std_com[compare_index,0],up_cyc_c[compare_index,0])
                            PCY = A_valuevel(mean_com[compare_index,1],std_com[compare_index,1],up_cyc_c[compare_index,1])

                            PRX = A_valuevel(mean_rhl[compare_index,0],std_rhl[compare_index,0],up_cyc_r[compare_index,0])
                            PRY = A_valuevel(mean_rhl[compare_index,1],std_rhl[compare_index,1],up_cyc_r[compare_index,1])
                            PRZ = A_valuevel(mean_rhl[compare_index,2],std_rhl[compare_index,2],up_cyc_r[compare_index,2])

                            PLX = A_valuevel(mean_lhl[compare_index,0],std_lhl[compare_index,0],up_cyc_l[compare_index,0])
                            PLY = A_valuevel(mean_lhl[compare_index,1],std_lhl[compare_index,1],up_cyc_l[compare_index,1])
                            PLZ = A_valuevel(mean_lhl[compare_index,2],std_lhl[compare_index,2],up_cyc_l[compare_index,2])

                            PCXV = A_valuevel(mean_comv[compare_index,0],std_comv[compare_index,0],up_cyc_cv[compare_index,0])
                            PCYV = A_valuevel(mean_comv[compare_index,1],std_comv[compare_index,1],up_cyc_cv[compare_index,1])

                            PRXV = A_valuevel(mean_rhlv[compare_index,0],std_rhlv[compare_index,0],up_cyc_rv[compare_index,0])
                            PRYV = A_valuevel(mean_rhlv[compare_index,1],std_rhlv[compare_index,1],up_cyc_rv[compare_index,1])
                            PRZV = A_valuevel(mean_rhlv[compare_index,2],std_rhlv[compare_index,2],up_cyc_rv[compare_index,2])

                            PLXV = A_valuevel(mean_lhlv[compare_index,0],std_lhlv[compare_index,0],up_cyc_lv[compare_index,0])
                            PLYV = A_valuevel(mean_lhlv[compare_index,1],std_lhlv[compare_index,1],up_cyc_lv[compare_index,1])
                            PLZV = A_valuevel(mean_lhlv[compare_index,2],std_lhlv[compare_index,2],up_cyc_lv[compare_index,2])
                            
                            Psum = PCX+PCY+PRX+PRZ+PLZ+PLZ +PRXV  +PRZV +PLZV+PLZV +PCXV +PCYV +PLY+PRY+PLYV+PRYV
                            psum_cont.append(Psum)
                            if Psum > 0.125:
                                pertflag = 1
                                print("perted")
                                pind = ((compare_index/100)*cyc_len)+(t-cyc_len)
                                pert = 2005 - timet0[0]
                                if trial == 18:
                                    pert = 1510
                                x = ((pert - indxrt[np.where(indxrt<pert)[0][-1]])/(indxrt[np.where(indxrt<pert)[0][-1]+1]-indxrt[np.where(indxrt<pert)[0][-1]]))*100
                                delay = (pind - pert)/100
                                X.append(x)
                                
                                if delay < 0 and abs(delay)<10:
                                    delay = abs(delay)

                                    
                                if t>pert and delay >-10:
                                    print(delay)
                                    avgd.append(delay)
                                    if early_flag == 0:
                                        good_det +=1 
                                    else:
                                        false_det +=1
                                        print("det after early")
                                    early_flag = 0
                                    compared_cyc.append(heel_strike_count-14)
                                    break
                                else: 
                                    early_flag = 1
                                    print(heel_strike_count)
                                 
                            else: 
                                pertflag = 0
                        if not pertflag:
                            #if not perturbed then update mean 
                            mean_com = np.mean((np.array(c_s)[-10:,:,:]),axis=0)
                            std_com = np.std((np.array(c_s)[-10:,:,:]),axis=0)
                            mean_rhl = np.mean((np.array(r_s)[-10:,:,:]),axis=0)
                            std_rhl = np.std((np.array(r_s)[-10:,:,:]),axis=0)
                            mean_lhl = np.mean((np.array(l_s)[-10:,:,:]),axis=0)
                            std_lhl = np.std((np.array(l_s)[-10:,:,:]),axis=0)
                            mean_comv = np.mean((np.array(c_s_v)[-10:,:,:]),axis=0)
                            std_comv = np.std((np.array(c_s_v)[-10:,:,:]),axis=0)
                            mean_rhlv = np.mean((np.array(r_s_v)[-10:,:,:]),axis=0)
                            std_rhlv = np.std((np.array(r_s_v)[-10:,:,:]),axis=0)
                            mean_lhlv = np.mean((np.array(l_s_v)[-10:,:,:]),axis=0)
                            std_lhlv = np.std((np.array(l_s_v)[-10:,:,:]),axis=0)
                        else:
                            break
                    #clear values to start next cycle 
                    cyc_c = []
                    cyc_l = []
                    cyc_r = []

                    cyc_cv = []
                    cyc_lv = []
                    cyc_rv = []

                tracking_active = True 
            if tracking_active: 
                #update cycle vector
                cyc_c.append(com_rel)
                cyc_l.append(lhl_rel)
                cyc_r.append(rhl_rel)

                cyc_cv.append(cur_v_c)
                cyc_lv.append(cur_v_l)
                cyc_rv.append(cur_v_r)
print("Average Delay(s)")
print(np.mean(np.array(avgd)))
print("Good")
print(good_det)
print("False")
print(false_det)
print("Pert occurance")
print(X)
print(compared_cyc)
total = np.sum(np.array(compared_cyc))
print(total)
 
# #11/14 updates mean after each cycle and not just comapred to pert

# import numpy as np
# from scipy.interpolate import interp1d
# import matplotlib.pyplot as plt
# import pickle
# #Load the dictionary from the pickle file
# with open('data_stream_AB03', 'rb') as f:
#     loaded_data = pickle.load(f)

# X = []
# false_det = 0
# no_det = 0

# pert = 2000

# # Access the variables
# time = loaded_data['time']
# COM_ALL = loaded_data['COM_ALL']
# RHL_ALL = loaded_data['RHL_ALL']
# LHL_ALL = loaded_data['LHL_ALL']
# VEL_COM_ALL = loaded_data['VEL_COM_ALL']
# VEL_RHL_ALL = loaded_data['VEL_RHL_ALL']
# VEL_LHL_ALL = loaded_data['VEL_LHL_ALL']
# indxr = loaded_data['indxr']
# psum_cont = []

# last_force = 0
# heel_strike = False
# frame = 0
# heel_strike_count = 0
# cont_c_rel = []
# cont_r_rel = []
# cont_l_rel = []
# #raw cont
# cont_c_raw = []
# cont_r_raw = []
# cont_l_raw = []
# #ony current cycle values
# cyc_c = []
# cyc_r = []
# cyc_l = []

# #raw current cycle
# cyc_c_raw = []
# cyc_r_raw = []
# cyc_l_raw = []
# cyc_len_list = []
# #split cycles, cont
# c_s = []
# r_s = []
# l_s = []
# #velocity over each cycle 
# c_s_v = []
# r_s_v =[]
# l_s_v = []
# #only this cycle
# cyc_cv = []
# cyc_rv = []
# cyc_lv = []



# getted_frames_cycle = []
# pert_trigger = False
# tracking_active = False
# weight_sum = False
# early_flag = 0
# def A_value(cpmean,cpstd,cpoint):
#     A = 0
#     lambdaa= 0
#     if (cpmean-2*cpstd)>cpoint:
#         #cpstd = max(cpstd,0.001)
#         A = abs(cpoint-(cpmean-2*cpstd)) 
#         lambdaa = (1/16)*(1/(2*abs(cpstd)+A))
#     if (cpmean+2*cpstd)<cpoint:
       
#         A = abs(cpoint-(cpmean+2*cpstd))
#         lambdaa = (1/16)*(1/(2*abs(cpstd)+A))
#     Pm = np.multiply(A,lambdaa)

#     return Pm

# def A_valuevel(cpmean,cpstd,cpoint):
#     A = 0
#     lambdaa= 0
#     if (cpmean-2*cpstd)>cpoint:
#         A = abs(cpoint-(cpmean-2*cpstd)) 
#         cpstd = cpstd/cpmean
#         lambdaa = (1/16)*(1/(2*abs(cpstd)+A))
#     if (cpmean+2*cpstd)<cpoint:
#         A = abs(cpoint-(cpmean+2*cpstd))
#         cpstd = cpstd/cpmean
#         lambdaa = (1/16)*(1/(2*abs(cpstd)+A))
#     Pm = np.multiply(A,lambdaa)

#     return Pm

# def resample(positions, num_points=100):
#     """
#     Resample foot position data to be evenly spaced between 0 and 100% of the gait cycle.

#     Parameters:
#     foot_positions (list or np.ndarray): A 2D list or array of shape (N, 3), where N is the number of original timestamps,
#                                          and each row represents the (x, y, z) position of the foot.
#     num_points (int): The number of points to interpolate to, representing the 0–100% gait cycle.

#     Returns:
#     np.ndarray: A 2D array of shape (num_points, 3) containing the resampled foot positions.
#     """
#     positions = np.array(positions)  # Ensure input is an array
#     num_original_points = positions.shape[0]
   
#     # Original timestamps (normalized from 0 to 100% based on the number of input points)
#     original_timestamps = np.linspace(0, 100, num_original_points)
   
#     # Target timestamps for interpolation (0% to 100% in the specified number of points)
#     target_timestamps = np.linspace(0, 100, num_points)
   
#     # Interpolate each axis (x, y, z) independently
#     x_interp = interp1d(original_timestamps, positions[:, 0])
#     y_interp = interp1d(original_timestamps, positions[:, 1])
#     z_interp = interp1d(original_timestamps, positions[:, 2])
   
#     # Generate resampled points for each axis at the target timestamps
#     resampled_x = x_interp(target_timestamps)
#     resampled_y = y_interp(target_timestamps)
#     resampled_z = z_interp(target_timestamps)
   
#     # Combine into a 2D array for the full resampled cycle
#     resampled_positions = np.vstack((resampled_x, resampled_y, resampled_z)).T
   
#     return resampled_positions

# trials = len(COM_ALL)
# avgd = []
# good_det = 0
# false_det = 0
# pertflag = 0
# for trial in range(trials):
#     #plt.plot(list(range(100)),np.swapaxes(np.array(r_s_v)[:,:,1],0,1),color="blue");plt.plot(mean_rhlv[:,1]+2*std_rhlv[:,1],color = 'black');plt.plot(mean_rhlv[:,1]-2*std_rhlv[:,1],color = 'black');plt.plot(mean_rhlv[:,1],color = 'black');plt.plot(list(range(100)),np.array(up_cyc_rv)[:,1],color="red");plt.plot(list(range(100)),np.array(r_s_v)[11,:,1],color="green")
#     heel_strike_count = 0
#     c_s = []
#     r_s = []
#     l_s = []
#     c_s_v = []
#     r_s_v = []
#     l_s_v = []
#     pind = 0
#     early_flag = 0
#     delay = 0
#     # cyc_c = []
#     # cyc_l = []
#     # cyc_r = []

#     # cyc_cv = []
#     # cyc_lv = []
#     # cyc_rv = []
#     timet0 = time[trial]
#     timet = timet0 - timet0[0]
#     indxrt = indxr[trial]
#     comt = COM_ALL[trial]
#     rhlt = RHL_ALL[trial]
#     lhlt = LHL_ALL[trial]
#     comvelt = VEL_COM_ALL[trial]
#     rhlvelt = VEL_RHL_ALL[trial]
#     lhlvelt = VEL_LHL_ALL[trial]
#     for ind,t in enumerate(timet):
#         if pert == 1:
#             break
#         com_rel = comt[ind,:]
#         rhl_rel = rhlt[ind,:]
#         lhl_rel = lhlt[ind,:]
#         if t > timet[0]:
#             cur_v_c = comvelt[ind,:]
#             cur_v_r = rhlvelt[ind,:]
#             cur_v_l = lhlvelt[ind,:]
#             if heel_strike_count >= len(indxrt):
#                     break
#             heel_strike = 1 if t == indxrt[heel_strike_count] else 0
#             if heel_strike:
#                 heel_strike_count +=1
#                 if heel_strike_count >= len(indxrt):
#                     break
#                 print("strike")
#                 if tracking_active:
#                     cyc_len,_ = np.shape(cyc_c)
#                     cyc_len_list.append(cyc_len)
#                     #mean and std first 10 cycles
#                     if heel_strike_count == 4:
#                         #COM
#                         up_cyc_c =resample(np.array(cyc_c))
#                         c_s.append(up_cyc_c)
#                         #plt.plot(list(range(100)),np.swapaxes(np.array(c_s)[:,:,0],0,1))
#                         mean_com = up_cyc_c
#                         #RIGHT
#                         up_cyc_r =resample(np.array(cyc_r))
#                         r_s.append(up_cyc_r)
#                         mean_rhl = up_cyc_r
#                         #LEFT
#                         up_cyc_l =resample(np.array(cyc_l))
#                         l_s.append(up_cyc_l)
#                         mean_lhl = up_cyc_l
#                         #append to list for making a continous mean 

#                         #Velocity 
#                         up_cyc_cv =resample(cyc_cv)
#                         #add to spi cycles vector for later mean
#                         c_s_v.append(up_cyc_cv)
#                         mean_comv = up_cyc_cv

#                         up_cyc_rv =resample(cyc_rv)
#                         r_s_v.append(up_cyc_rv)
#                         mean_rhlv = up_cyc_rv 

#                         up_cyc_lv =resample(cyc_lv)
#                         l_s_v.append(up_cyc_lv)
#                         mean_lhlv = up_cyc_lv 
                        
#                     if ( heel_strike_count >4  and heel_strike_count<=14):
#                         up_cyc_c = resample(np.array(cyc_c))
#                         c_s.append(up_cyc_c)
#                         mean_com = np.mean((np.array(c_s)),axis=0)
#                         std_com = np.std((np.array(c_s)),axis=0)

#                         up_cyc_r = resample(np.array(cyc_r))
#                         r_s.append(up_cyc_r)
#                         mean_rhl = np.mean((np.array(r_s)),axis=0)
#                         std_rhl = np.std((np.array(r_s)),axis=0)

#                         up_cyc_l = resample(np.array(cyc_l))
#                         l_s.append(up_cyc_l)
#                         mean_lhl = np.mean((np.array(l_s)),axis=0)
#                         std_lhl = np.std((np.array(l_s)),axis=0)
                        
#                         #Velocity 
#                         up_cyc_cv =resample(cyc_cv)
#                         #add to spi cycles vector for later mean
#                         c_s_v.append(up_cyc_cv)
#                         mean_comv = np.mean((np.array(c_s_v)),axis=0)
#                         std_comv = np.std((np.array(c_s_v)),axis=0)

#                         up_cyc_rv =resample(cyc_rv)
#                         r_s_v.append(up_cyc_rv)
#                         mean_rhlv = np.mean((np.array(r_s_v)),axis=0)
#                         std_rhlv = np.std((np.array(r_s_v)),axis=0)

#                         up_cyc_lv =resample(cyc_lv)
#                         l_s_v.append(up_cyc_lv)
#                         mean_lhlv = np.mean((np.array(l_s_v)),axis=0)
#                         std_lhlv = np.std((np.array(l_s_v)),axis=0)

#                     if heel_strike_count>14:
#                         if early_flag ==1:
#                             break
#                         weight_sum = True
#                         up_cyc_c = resample(np.array(cyc_c))
#                         c_s.append(up_cyc_c)
                        
#                         up_cyc_r = resample(np.array(cyc_r))
#                         r_s.append(up_cyc_r)
                        
#                         up_cyc_l = resample(np.array(cyc_l))
#                         l_s.append(up_cyc_l)
#                         #Velocity 
#                         up_cyc_cv =resample(cyc_cv)
#                         #add to spi cycles vector for later mean
#                         c_s_v.append(up_cyc_cv)
#                         up_cyc_rv =resample(cyc_rv)
#                         r_s_v.append(up_cyc_rv)
                        
#                         up_cyc_lv =resample(cyc_lv)
#                         l_s_v.append(up_cyc_lv)
#                         psum_cont = []
#                         for compare_index in range (100):
#                             PCX = A_valuevel(mean_com[compare_index,0],std_com[compare_index,0],up_cyc_c[compare_index,0])
#                             PCY = A_valuevel(mean_com[compare_index,1],std_com[compare_index,1],up_cyc_c[compare_index,1])

#                             PRX = A_valuevel(mean_rhl[compare_index,0],std_rhl[compare_index,0],up_cyc_r[compare_index,0])
#                             PRY = A_valuevel(mean_rhl[compare_index,1],std_rhl[compare_index,1],up_cyc_r[compare_index,1])
#                             PRZ = A_valuevel(mean_rhl[compare_index,2],std_rhl[compare_index,2],up_cyc_r[compare_index,2])

#                             PLX = A_valuevel(mean_lhl[compare_index,0],std_lhl[compare_index,0],up_cyc_l[compare_index,0])
#                             PLY = A_valuevel(mean_lhl[compare_index,1],std_lhl[compare_index,1],up_cyc_l[compare_index,1])
#                             PLZ = A_valuevel(mean_lhl[compare_index,2],std_lhl[compare_index,2],up_cyc_l[compare_index,2])

#                             PCXV = A_valuevel(mean_comv[compare_index,0],std_comv[compare_index,0],up_cyc_cv[compare_index,0])
#                             PCYV = A_valuevel(mean_comv[compare_index,1],std_comv[compare_index,1],up_cyc_cv[compare_index,1])

#                             PRXV = A_valuevel(mean_rhlv[compare_index,0],std_rhlv[compare_index,0],up_cyc_rv[compare_index,0])
#                             PRYV = A_valuevel(mean_rhlv[compare_index,1],std_rhlv[compare_index,1],up_cyc_rv[compare_index,1])
#                             PRZV = A_valuevel(mean_rhlv[compare_index,2],std_rhlv[compare_index,2],up_cyc_rv[compare_index,2])

#                             PLXV = A_valuevel(mean_lhlv[compare_index,0],std_lhlv[compare_index,0],up_cyc_lv[compare_index,0])
#                             PLYV = A_valuevel(mean_lhlv[compare_index,1],std_lhlv[compare_index,1],up_cyc_lv[compare_index,1])
#                             PLZV = A_valuevel(mean_lhlv[compare_index,2],std_lhlv[compare_index,2],up_cyc_lv[compare_index,2])
                            
#                             Psum = PCX+PCY+PRX+PRZ+PLZ+PLZ +PRXV  +PRZV +PLZV+PLZV +PCXV +PCYV +PLY+PRY+PLYV+PRYV
#                             psum_cont.append(Psum)
#                             if Psum > 0.12:
#                                 pertflag = 1
#                                 print("perted")
#                                 pind = ((compare_index/100)*cyc_len)+(t-cyc_len)
#                                 pert =  2005 - timet0[0]
#                                 if trial ==18:
#                                     pert = 1510
#                                 x = ((pert - indxrt[np.where(indxrt<pert)[0][-1]])/(indxrt[np.where(indxrt<pert)[0][-1]+1]-indxrt[np.where(indxrt<pert)[0][-1]]))*100
#                                 delay = (pind - pert)/100
#                                 X.append(x)
                                
#                                 if delay < 0 and abs(delay)<10:
#                                     delay = abs(delay)
                                    
                                    
#                                 if t>pert and delay >-10:
#                                     print(delay)
#                                     avgd.append(delay)
#                                     good_det +=1 
#                                     early_flag = 0
#                                     break
#                                 else: 
#                                     early_flag = 1
#                                     false_det +=1
#                                     print(heel_strike_count)
#                                     break
#                             else: 
#                                 pertflag = 0
#                         if not pertflag:
#                             #if not perturbed then update mean 
#                             mean_com = np.mean((np.array(c_s)[-10:,:,:]),axis=0)
#                             std_com = np.std((np.array(c_s)[-10:,:,:]),axis=0)
#                             mean_rhl = np.mean((np.array(r_s)[-10:,:,:]),axis=0)
#                             std_rhl = np.std((np.array(r_s)[-10:,:,:]),axis=0)
#                             mean_lhl = np.mean((np.array(l_s)[-10:,:,:]),axis=0)
#                             std_lhl = np.std((np.array(l_s)[-10:,:,:]),axis=0)
#                             mean_comv = np.mean((np.array(c_s_v)[-10:,:,:]),axis=0)
#                             std_comv = np.std((np.array(c_s_v)[-10:,:,:]),axis=0)
#                             mean_rhlv = np.mean((np.array(r_s_v)[-10:,:,:]),axis=0)
#                             std_rhlv = np.std((np.array(r_s_v)[-10:,:,:]),axis=0)
#                             mean_lhlv = np.mean((np.array(l_s_v)[-10:,:,:]),axis=0)
#                             std_lhlv = np.std((np.array(l_s_v)[-10:,:,:]),axis=0)
#                         else:
#                             break
#                     #clear values to start next cycle 
#                     cyc_c = []
#                     cyc_l = []
#                     cyc_r = []

#                     cyc_cv = []
#                     cyc_lv = []
#                     cyc_rv = []

#                 tracking_active = True 
#             if tracking_active: 
#                 #update cycle vector
#                 cyc_c.append(com_rel)
#                 cyc_l.append(lhl_rel)
#                 cyc_r.append(rhl_rel)

#                 cyc_cv.append(cur_v_c)
#                 cyc_lv.append(cur_v_l)
#                 cyc_rv.append(cur_v_r)
# print("Average Delay (s)")
# print(np.mean(np.array(avgd)))
# print("Good")
# print(good_det)
# print("False")
# print(false_det)
# print("Pert occurance")
# print(X)
                
