#SORTING DATA BASED ON % GAIT CYCLE 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import os 
from scipy.signal import decimate
from mpl_toolkits.mplot3d import Axes3D




# #define all functions ------------------------------------------------------------------------
# def interp_nan(data):
#     """
#     Interpolates NaN values in each column of a 2D array separately.

#     Parameters:
#         data (np.ndarray): 2D NumPy array where NaNs will be interpolated in each column.

#     Returns:
#         np.ndarray: A 2D array with NaNs replaced by interpolated values.
#     """
    # # Make a copy to avoid modifying the original data
    # interpolated_data = np.array(data, copy=True)
    
    # # Interpolate NaN values in each column
    # for col in range(interpolated_data.shape[1]):
    #     # Find indices where values are not NaN
    #     nans = np.isnan(interpolated_data[:, col])
    #     not_nans = ~nans
        
    #     # If the entire column is NaN, skip interpolation
    #     if not np.any(not_nans):
    #         continue
        
    #     # Perform linear interpolation over NaNs in the current column
    #     interpolated_data[nans, col] = np.interp(
    #         np.flatnonzero(nans), 
    #         np.flatnonzero(not_nans), 
    #         interpolated_data[not_nans, col]
    #     )
    
    # return interpolated_data


def interp_nans(data):
    """
    Interpolates NaN values along each column of a 2D array.
    
    Parameters:
    data (numpy.ndarray): The input 2D data array containing NaNs.
    
    Returns:
    numpy.ndarray: The data array with NaNs linearly interpolated along each column.
    """
    if data.ndim != 2:
        raise ValueError("Input data must be a 2D array.")
    
    # Convert the data into a pandas DataFrame
    df = pd.DataFrame(data)
    
    # Apply interpolation along each column
    interpolated_df = df.interpolate(method='linear', axis=0, limit_direction='both')
    
    # Convert back to a numpy array
    return interpolated_df.to_numpy()

#LOAD ALL DATA FOR ONE PERT TYPE FOR ONE PARTICIPANT FROM THE SPREADSHEET
def dataloader(subject, trial):
     #load filepath and data for markers
    marker_dir = os.path.join(r"C:\Users\maria\OneDrive\Desktop\PhD_re\ct302024", subject,'Marker',trial)
    
    marker_table = pd.read_csv(marker_dir)
    
    grf_dir = os.path.join(r"C:\Users\maria\OneDrive\Desktop\PhD_re\ct302024", subject, 'GRF',trial)

    grf_table = pd.read_csv(grf_dir)
    #convert to np array
    mtt = marker_table.values
    grf = grf_table.values
    #extract time data
    time_vec0 = mtt[:, 0]
    
    mt = interp_nans(mtt)
        
    #Pull correct columns
    pelx = np.column_stack([mt[:, 14], mt[:, 17], mt[:, 20], mt[:, 23]])
    COMx = -np.mean(pelx, axis=1)

    pely = np.column_stack([mt[:, 15], mt[:,18], mt[:, 21], mt[:, 24]])
    COMy = np.mean(pely, axis=1)

    pelz = np.column_stack([mt[:, 16], mt[:, 19], mt[:, 22], mt[:, 25]])
    COMz = np.mean(pelz, axis=1)

    COM = np.column_stack([COMy, COMx, COMz])

    RHLx = mt[:, 8]
    RHLy = mt[:, 9]
    RHLz = mt[:, 10]
    
    RHL = np.column_stack([RHLy, RHLx, RHLz])

    LHLx = mt[:, 2] #23
    LHLy = mt[:, 3] #24
    LHLz = mt[:, 4] #25
    
    LHL = np.column_stack([LHLy, LHLx, LHLz])

    Rtoex = mt[:, 11]
    Rtoey = mt[:, 12]
    Rtoez = mt[:, 13]

    Rtoe = np.column_stack([Rtoey, Rtoex, Rtoez])
    
    Ltoex = mt[:, 5] #20
    Ltoey = mt[:, 6] #21
    Ltoez = mt[:, 7] #22

    Ltoe = np.column_stack([Ltoey, Ltoex, Ltoez])

    RGRF = abs(grf[:, 4])
    LGRF = abs(grf[:, 13])
    
    #downsample grf to match marker data frequency/size
    downsample_factor = 10
    RGRF = RGRF[:-1]
    LGRF = LGRF[:-1]
    
    # Downsample the data using decimate
    RGRF = decimate(RGRF, downsample_factor, ftype='iir')
    LGRF = decimate(LGRF, downsample_factor, ftype='iir')
    
    
    #size data points, 3 (xyz), 5
    trial_data = np.concatenate((COM[:,:,np.newaxis],RHL[:,:,np.newaxis],LHL[:,:,np.newaxis],Rtoe[:,:,np.newaxis],Ltoe[:,:,np.newaxis]),axis=2)
    return COM/1000, RHL/1000, LHL/1000, Rtoe/1000, Ltoe/1000, RGRF, LGRF,time_vec0, trial_data/1000

#DETERMINE THE LOCATIONS OF THE GAIT CYCLES USING GRF DATA
def grf_strike(GRF):
    threshold = 70
    bi_GRF = GRF.copy()
    
    # Binarize the ground reaction forces
    bi_GRF[bi_GRF < threshold] = 0
    bi_GRF[bi_GRF >= threshold] = 1  # Changed to >= to catch the exact threshold
    
    # Calculate the difference
    diff_bi = np.diff(bi_GRF)
    
    # Find indices where a heel strike occurs (transition from 0 to 1)
    indicies = np.where(diff_bi == 1)[0] 
    
    # Initialize list to store valid heel strike indices
    valid_strikes = []
    last_strike_index = -80  # Initialize to a value that allows the first strike
    
    # Filter the indices based on the distance condition
    for index in indicies:
        if (index - last_strike_index) >= 80:
            valid_strikes.append(index)
            last_strike_index = index  # Update the last detected index
    
    return np.array(valid_strikes)

#DETERMINE THE LOCATIONS OF PERT
def pert_indexx(time_vec0):
    pert_ind = np.where(time_vec0==2000)
    return pert_ind[0][0]


#SPLIT THE ENTIRE DATASET INTO CYCLES BASED ON THE HEEL STRIKE LOCATIONS AND UPSAMPLE RELATIVE TO GAIT PHASE

#FIND THE VELOCITY IN X,Y,Z ACROSS THE DATA AND UPSAMPLED TO GAIT PHASE  
def velocity(data,time):
       #Calculate Velocity
        velocity = np.zeros((len(data[:,0])-1, 3))
        velocity[:,0] = np.divide(np.diff(data[:,0]),np.diff(time))
        velocity[:,1] = np.divide(np.diff(data[:,1]),np.diff(time))
        velocity[:,2] = np.divide(np.diff(data[:,2]),np.diff(time))

        return velocity


#NORMALIZE ALL DATA RELATIVE TO THE INITIAL VALUE OF THE COM, RIGHT AND LEFT HEEL 
def norm_rel(trial_data):
    #distance from right foor to com
    RHL_rel = trial_data[:,:,1]-trial_data[:,:,0]
    #distance from left foot to com
    LHL_rel = trial_data[:,:,2]-trial_data[:,:,0]
    MIDPOINT_FEET = ((trial_data[:,:,1]+trial_data[:,:,2]+trial_data[:,:,3]+trial_data[:,:,4]))/4
    COM_rel = trial_data[:,:,0] - MIDPOINT_FEET
    return COM_rel, RHL_rel, LHL_rel

#INITIALIZE ARRAYS TO STORE ALL THE OUTPUTS 
COM_ALL = []
RHL_ALL = []
LHL_ALL = []

VEL_COM_ALL =[]
VEL_RHL_ALL = []
VEL_LHL_ALL = []
time = []
X = []

INDXR = []
filesab04 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
#filesab04edit = [1,2,3,4,5,6,9,11,12,13,16,17]
for i in filesab04:
    #generate the data for all three trials of subject 
    if i <=10:
        filename1 = f"OT_{i}.csv"
    else:
        filename1 = f"OS_{i-10}.csv"
    
    # filename1 = f"T_0{i}.csv"
    COM, RHL, LHL, Rtoe, Ltoe, RGRF, LGRF, time_vec0,trial_data = dataloader(subject = "AB04",trial = filename1)
    
    #find index of heel strikes and output a matrix of all the data of the trial
    indxR = grf_strike(RGRF)
    #plt.plot(time_vec0,RGRF[:-1]);plt.plot(time_vec0[indxR],RGRF[indxR],'o')

    cycle_length = int(np.round(np.mean(np.abs(np.diff(indxR[4:12])))))
    
    time.append(time_vec0)
    num_r = len(indxR)
    num_r_af_14 = num_r-12
    rr = np.array(list(range(12, num_r)))
    for r in rr:
        if indxR[r]-indxR[r-1] < cycle_length - 15 or indxR[r]-indxR[r-1] > cycle_length + 15:
            indxR[r]=indxR[r-1]+cycle_length
    if indxR[-1]> len(COM):
        indxR = np.delete(indxR,-1)
    COM_norm, RHL_rel,LHL_rel= norm_rel(trial_data)

    vel_com = velocity(trial_data[:,:,0],time_vec0)
    vel_rhl = velocity(trial_data[:,:,1],time_vec0)
    vel_lhl = velocity(trial_data[:,:,2],time_vec0)
 
    VEL_COM_ALL.append(vel_com)
    VEL_RHL_ALL.append(vel_rhl)
    VEL_LHL_ALL.append(vel_lhl)
    
    #append to all previous ones
    COM_ALL.append(COM_norm)
    RHL_ALL.append(RHL_rel)
    LHL_ALL.append(LHL_rel)
    INDXR.append(indxR)
#AFTER LOOPING 
    
# Save variables to .npy files
import pickle

# Create a dictionary to store all the variables without Data_setup prefix
data_dict = {
    'time': time,
    'COM_ALL': COM_ALL ,
    'RHL_ALL': RHL_ALL,
    'LHL_ALL': LHL_ALL,
    'VEL_COM_ALL': VEL_COM_ALL,
    'VEL_RHL_ALL': VEL_RHL_ALL,
    'VEL_LHL_ALL': VEL_LHL_ALL,
    'indxr' :INDXR


}

# Save the dictionary to a pickle file
with open('data_stream_AB04', 'wb') as f:
    pickle.dump(data_dict, f)
