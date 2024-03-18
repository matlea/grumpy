"""

Methods for loading, manipulating, and displaying Prodigy data.

This is a major re-write of grump.py version 23.05.13, mainly regarding
the data dicts.
grump.py, up to version 23.05.13, was previously called GrumpyTwo.py but 
then renamed grumpy.py (from 23.03.13), replacing the original grumpy.py
that might still be floating around. 

Major re-write from 23.05.14. Not fully tested yet. Mostly tested. Pretty much fully tested.

In general: Intensities are stored in key 'int'. See below for specifics regarding spin data.
            Axes are: 'x' for 1d data, 'x' and 'y' for 2d data, and 'x', 'y', and 'z' for 3d data.
            Note that for e.g. line scans where both deflectors are used there is an 'y2' axis as well.
            Axis types and labels are found in the 'Meta' dict, as well as the intensity label.
            Key 'Measurement_type' describes the type of data, while key 'Type' is updated in case of 
            data manipulation.
            Experimental settings are found key 'Experiment', parameters and their values are found in
            keys 'Parameters' and 'Parameter_values', the raw data is found in key Data.
            About spin data: (description pending)

Measurement types:
    ARPES                                [y,x] = [angle, energy]
    XPS                                    [x] = [energy]
    spin_edc                               [x] = [energy]
    target_scattering_spectrum             [x] = [scattering energy]
    fermi_map                          [z,y,x] = [delf_x, angle, energy]
    scattering_intensity_scan              [x] = [time]
    deflector_scattering_intensity_scan  [y,x] = [deflection, energy]
    deflector_total_intensity_scan       [y,x] = [deflection, energy]
    spin_deflector_1a (1)                [y,x] = [deflection, energy]   
                                           [x] = [deflection]           for fixed energy
    spin_deflector_1b (1)                [y,x] = [deflection, energy]
    manipulator_line_scan                  [x] = [poisition]            for fixed energy

(1) Line scan with one deflector, (2) Line scan using both deflectors. 

Main methods:
    Load()
    Plot()
    SubArray()          Does not work for deflector spin data yet.
    Compact()           Does not work with spin data yet.
    Profile()           (combines SubArray() and Compact()).
    QuickSpin()
    SpinEDC()
    QuickSpinMDC()
    Polarization()      Not checked after re-write.
    MergeSpinEDC()      Merge two spin edc data sets into one (merged_data = MergeSpinEDC(data1,data2))
    RemoveSpinEDC()     Remove a specific curve from a spin edc dict.

    Explore()           For 2d and 3d data. Not quite ready but it works. Need to add option to rotate the images.
    ShiftAxis()
    PlotFermiMap()      (Is called by Plot() when needed)
    Fit()               Work in progress. Ready for 1d data with profiles gauss and gauss2.
    Smooth()            Work in progress. Ready for 1d and 2d data.
    PlotFitRes()        Plots fit results (also called by Plot() if needed).
    info()
    Info()
    SaveSpinEDC2File()  Save one or several spin edc data sets to text file.
    SaveARPES2File()
    SaveFermiMap2File()
    SaveTargetScatteringSpectrum2File()
    SaveXPS2File()
    Save()              Generic method for all the Save...() methods above.

    AppendFEmaps()      An unsophiticated method to append two Fermi maps.

    xps()               Not a method but a class, containing binding energies for elements.

    SortP()             (Loader for general data not covered by the types listed above. A work in progress. Slow work, slow progress.)

    ExportSpinEDC()     Save spin edc as a 3-column file (E, NegPol off, NegPol on)
    
    MakeArithmetic()    Makes an object out of a Grumpy dict which can handle arithmetic operations. Applies to non-spin data.

More methods (not applicapble on grumpy dicts):
    fitGauss()          Fit a gaussian.
    fitGauss2()         Fit two gaussians.
    medianFilter()      Run a median filter over an intensity curve. medianFilter(Y, size, mode).
    extendFit()         Applicable to fit result dicts.
    removeSpikes()      Remove spikes. removeSpikes(Y, threshold). Legacy name: filterOutliers()
    fwhm()              Simple algorithm to find fwhm.

Notes:
    1.  Methods starting with upper case accepts grumpy dicts, while methods starting with lower case
        accepts either arrays or fit dicts. Exception: info() accepts grumpy dicts.
    2.  To list a method's recognized kwargs, pass kw = True.

--- 

Version history (from 23.05.15 and onwards):

Version 24.03.18b   Added more save-to-file methods (SaveXPS2File() and SaveTargetScatteringSpectrum2File()) and included them in Save().
                    Changed all Fore.BLACK to Fore.RESET.

Version 24.03.17    Added SaveFermiMap2File() and included it in Save().

Version 24.03.16    Added SaveSpinEDC2File() (it seems like I already had made a ExportSpinEDC() method but this new one is different 
                    as it saves all the individual spin edc's, while the previos one only saves the average intensities.)
                    Added RemoveSpinEDC().
                    Added SaveARPES2File().
                    Started Save(). Save() is a method using other save-to-text methods. So far: SaveARPES2File() and SaveSpinEDC2File()

Version 24.03.15    Added new kwargs to QuickSpin. Normalize the edc's to each other before calculating the asymmetry. Kwargs are norm = True,
                    and norm_start_point and norm_points (e.g. 2 and 5 will use points 2 to 7 for normalization).

Version 24.01.24    Added 'Analyer CCD' as a recognized alternative to 'PhoibosCCD'. 'Analyzer CCD' was added to our Prodigy installation
                    during Robert's visit in January (2024).
                    NOTE: Add 'Analyzer Intensity' and 'Analyzer Spin' as well a.s.a.p. 

Version 23.12.16    Updated SpinEDC().
                    Updated MassageSpinEDC() (added some time ago but failed to update history then)

Version 23.12.16    Updated ExportSpinEDC(). This method was added a week ago but I forgot to update the history.
                    Added SpinEDC(). A method to manipulated spin edc data (normalize, exclude, etc) and create a new data dict.

Version 23.12.04    Added MergeSpinEDC() to merge spin edc:s with the same settings.

Version 23.11.20    Added a kwarg in Load() (xps_sum = False) that prevents data in HM2 mode (or similar) to be
                    compacted to one axis. Pass xps_sum = True to override this (i.e. True compacts).

Version 23.11.15    Added a method (MakeArithmetic) to make an objectc out of a Grumpy dict that can handle arithmetic operations.
                    Does not apply to spin data (yet).

Version 23.10.27    Added to gitlab

Version 23.09.27    Added AppendFEmaps() to append two femi maps. It comes with restrictions...

Version 23.08.06    Moved fwhm() from myRay.py to here.

Version 23.06.05    Added helper method _kwarg_checker() to be used by all methods using kwargs.
                    Added to _kwarg_checker(): if a kwarg kw = True is passed then all recognized
                    kwargs are listed.
                    Added Smooth().
                    Updated removeSpikes() to accept 2d data.

Version 23.06.04    Removed fitGauss_lb() & _gauss_lb() and changed fitGauss() and _gauss() to include
                    a linear background.
                    Added fitGauss2() and _gauss2() for double gaussians.
                    Added PlotFitRes(). Modified Plot() to call it when needed.

Version 23.06.02    Added class xps. Contains binding elergies for elements.

Version 23.06.01    Added method fitGauss_lb() to Fit().
                    Renamed filterOutliers() to removeSpikes().
                    Added extendFit(). 

Version 23.05.31    Added Fit(). Works for 1d and gauss. Add more.

Version 23.05.27    Added (re-added) method Profile(). It combines SubArray() and Compact().

Version 23.05.27    Buggfixes. Tests of method Polarization().

Version 23.05.23    General touch-ups.    
    
Version 23.05.21    Added Explore() for 2d and 3d data viewing. Not completely ready. There should be an option
                    to rotate images for 3d data (e.g. in case of Fermi maps).

Version 23.05.17    Improvements and buggfixes.
                    To do:  Fix SubArray() for deflector spin data.
                            Check Polarization().

Version 23.05.16    Finding and fixing buggs.

Version 23.05.15    Finished most of the re-writing of grumpy.py. The data format has been drastically simplified.
                    Last version with the old format is v.23.05.13.
                    Polarization() re-write is finished but the method is not tested.

"""

__version__ = "24.03.18"
__author__  = "Mats Leandersson"


import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit

from colorama import Fore
from sys import getsizeof
from copy import deepcopy

try: 
    import ipywidgets as ipw
    from IPython.display import display
except: 
    print(Fore.RED + '(\nCould not import the ipywidget module and/or display from IPython.display.') 
    print('Interactive plots will not work.\n\n' + Fore.RESET)

# ==============================================================================================
print(Fore.LIGHTWHITE_EX + f"grumpy, version {__version__}")
print("  SAL X = ShiftX, SALY = -ShiftY, Asymmetry: negative polarity off - on.")
print("  CCD: PhoibosCCD and AnalyzerCCD." + Fore.RESET)
# ==============================================================================================

# Globals

SHERMAN = 0.29

defax_energy = "energy"
defax_angle = 'angle'
defax_xdefl = 'deflector_x'
defax_ydefl = 'deflector_y'
defax_time = 'time'
defax_position = 'position'
defax_scattering_energy = 'scattering_energy'

# ==============================================================================================

CCD_ANALYZERS = ["PhoibosCCD", "AnalyzerCCD"]



# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================


def MeasurementType(data, shup = False, **kwargs):
    if not type(data) is dict:
        print(Fore.RED + 'Argument data must be a grumpy dict.' + Fore.RESET); return
    
    global CCD_ANALYZERS
    
    Experiment = data.get('Experiment', {})
    
    Measurement_type = ''
    #global Measurement_types
    
    ANALYZER = Experiment.get('Analyzer', '')
    NON_ENERGY_ORDINATE = np.array(data.get('Non_Energy_Ordinate', np.array([])))
    PARAMETERS = data.get('Parameters', [])
    

    # ARPES or XPS
    if len(PARAMETERS) == 0 and ANALYZER in CCD_ANALYZERS and len(NON_ENERGY_ORDINATE) > 0:
        if abs(NON_ENERGY_ORDINATE.min()) < 1 and abs(NON_ENERGY_ORDINATE.max()) < 1:
            Measurement_type = 'XPS'
        else:
            Measurement_type = 'ARPES'
    
    # Spin EDC (simple)
    elif len(PARAMETERS) == 1 and 'NegativePolarity' in PARAMETERS and ANALYZER == 'PhoibosSpin':
        Measurement_type = 'spin_edc'
    
    # Target scattering spectrum
    elif len(PARAMETERS) == 5 and "ScatteringEnergy [V]" in PARAMETERS and ANALYZER == 'PhoibosSpin':
        Measurement_type = 'target_scattering_spectrum'
    
    # Fermi map
    elif len(PARAMETERS) == 1 and ("SAL X [deg]" in PARAMETERS or "ShiftX [a.u.]" in PARAMETERS) and ANALYZER in CCD_ANALYZERS:
        Measurement_type = 'fermi_map'
    
    # Deflector scattering scan
    elif len(PARAMETERS) == 1 and ("SAL X [deg]" in PARAMETERS or "ShiftX [a.u.]" in PARAMETERS or
                                 "SAL Y [deg]" in PARAMETERS or "ShiftY [a.u.]" in PARAMETERS) and ANALYZER == 'PhoibosSpin':
        Measurement_type = 'deflector_scattering_intensity_scan'

    # Deflector Phoibos intensity scan
    elif len(PARAMETERS) == 1 and ("SAL X [deg]" in PARAMETERS or "ShiftX [a.u.]" in PARAMETERS or
                                 "SAL Y [deg]" in PARAMETERS or "ShiftY [a.u.]" in PARAMETERS) and ANALYZER == 'PhoibosIntensity':
        Measurement_type = 'deflector_total_intensity_scan'
    
    # Scattering scan
    elif len(PARAMETERS) == 0 and ANALYZER == 'PhoibosSpin':
        Column_labels = data.get('Experiment', {}).get('Column_labels', ['', ''])
        if Column_labels[0] == 'index':
            Measurement_type = 'scattering_intensity_scan'
    
    # Deflector spin
    elif ANALYZER == 'PhoibosSpin' and 'NegativePolarity' in PARAMETERS and len(PARAMETERS) > 1:
        if len(PARAMETERS) == 2:
            if ("SAL X [deg]" in PARAMETERS or "ShiftX [a.u.]" in PARAMETERS) ^ ("SAL Y [deg]" in PARAMETERS or "ShiftY [a.u.]" in PARAMETERS):
                Measurement_type = 'spin_deflector_1a'
        elif len(PARAMETERS) == 3:
            if ("SAL X [deg]" in PARAMETERS or "ShiftX [a.u.]" in PARAMETERS) and ("SAL Y [deg]" in PARAMETERS or "ShiftY [a.u.]" in PARAMETERS):
                Measurement_type = 'spin_deflector_1b'
                # here it is not clear if it is a 1d or 2d measurement. must address that later!
    
    # Omniax line scans
    elif len(PARAMETERS) == 1 and 'Position [mm]' in PARAMETERS:
        # This will have to be changed once we start doing more sophisticated manipulator stuff
        Measurement_type = 'manipulator_line_scan'


    if Measurement_type == '' and not shup:
        print('Could not identify the type of measurement.')
    
    return Measurement_type


# ==============================================================================================



def Load(file_name = '', path = '', shup = True, **kwargs):
    """
    Load exported data from Prodigy.
    Identifies common measurement setups.
    """
    recognized_kwargs = ['remove_spikes', 'filter_outliers', 'threshold', "sum_xps"]
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    DC = {}             # The main dict
    File = {}           # Sub dict, add to dict DC
    Experiment = {}     # Sub dict, add to dict DC
    Parameters = []     # Array, add to dict Experiment

    try: nshup = not shup
    except: nshup = False

    if not type(file_name) is str or not type(path) is str:
        print(Fore.RED + "Args file_name and path must be strings." + Fore.RESET)
        return {}
    if file_name == '':
        print(Fore.RED + "Arg file_name can not be an empty string." + Fore.RESET)
        return {}
    if not file_name.lower().endswith('xy'): file_name = f"{file_name}.xy"
    File.update({'file_name': file_name})
    File.update({'file_extension': file_name.split('.')[-1]})
    if path == '':
        File.update({'file_path': ''})
        File.update({'file': file_name})
    else:
        if path.endswith('/'): path = path[:-2]
        File.update({'file_path': path})
        File.update({'file': path + '/' + file_name})
    try:
        f = open(File.get('file'), 'r')
        f.close()
        File.update({'file_exists': True})
    except:
        File.update({'file_exists': False})
        print(Fore.RED + 'Could not find/open the file.' + Fore.RESET)
        return {}

    # ============= READ HEADER

    with open(File.get('file'), 'rb') as f:
        ReadHeader = True
        DataRowStart = 0
        
        #if nshup: print('reading header...')
        for i, row_ in enumerate(f):
            row = row_.decode(errors='ignore')
            row = row.replace('\r', '').replace('\n', '')

            if ReadHeader:
                if row.startswith('# Created by'):
                    Experiment.update({'Version': row.split(",")[1].replace(' ','').strip('Version')})
                if row.startswith('#   Energy Axis'):
                    Experiment.update({'Energy_Axis': row.split(":")[1].replace(' ','')})
                if row.startswith('#   Count Rate'):
                    Experiment.update({'Count_Rate': row.split(": ")[1].replace(' ','').replace('per','/')})
                if row.startswith('# Spectrum ID:'):
                    Experiment.update({'Spectrum_ID': int(row.split(":")[1].replace(' ',''))})
                if row.startswith("# Analysis Method"):   
                    Experiment.update({'Analysis_Method': row.split(":")[1].replace(' ','')})
                if row.startswith("# Analyzer:"):   
                    Experiment.update({'Analyzer': row.split(":")[1].replace(' ','')})
                if row.startswith("# Analyzer Lens:"):   
                    Experiment.update({'Lens_Mode': row.split(":")[1].replace(' ','')})
                if row.startswith("# Scan Mode:"):                                          
                    Experiment.update({'Scan_Mode': row.split(":")[1].replace(' ','')})
                if row.startswith("# Curves/Scan:"): 
                    Experiment.update({'Curves_Per_Scan': int(row.split(":")[1])})
                if row.startswith("# Values/Curve:"): 
                    Experiment.update({'Values_Per_Curve': int(row.split(":")[1])})
                if row.startswith("# Dwell Time:"): 
                    Experiment.update({'Dwell_Time': float(row.split(":")[1])})
                if row.startswith("# Excitation Energy:"):
                    Experiment.update({'Excitation_Energy': float(row.split(":")[1].replace(' ',''))})
                if row.startswith("# Kinetic Energy:"):   
                    Experiment.update({'Kinetic_Energy': float(row.split(":")[1])})
                    Experiment.update({'Ek': Experiment.get('Kinetic_Energy')})
                if row.startswith("# Pass Energy"): 
                    Experiment.update({'Pass_Energy': float(row.split(":")[1])})
                    Experiment.update({'Ep': Experiment.get('Pass_Energy')})
                if row.startswith("# OrdinateRange"):
                    tmp = row.split(":")[1].strip(' ').strip('[').strip(']').split(',')
                    Experiment.update({'Ordinate_Range': [float(tmp[0]), float(tmp[1])]})
                    Experiment.update({'Ordinate_Range_Start': float(tmp[0])})
                    Experiment.update({'Ordinate_Range_End': float(tmp[1])})
                if row.startswith("# Parameter:") and not "Step" in row:
                    par = row.split(":")[1].split('=')[0].replace('" ', '').replace(' "', '')
                    Parameters.append(par)
                if row.startswith("Number of Scans:"):
                    Experiment.update({'Number_of_scans': int(row.split(':')[1])})
                if row == '# Cycle: 0':
                    DataRowStart = i
                if row.startswith('# ColumnLabels'):
                    try:
                        column_labels = row.split(':')[1][17:].split(' ')
                    except:
                        column_labels = ['?', '?']
                    ReadHeader = False
            else:
                f.close()
                break
        Parameters = list(np.unique(Parameters))
        Experiment.update({'Parameters': Parameters})
        Experiment.update({'Column_labels': column_labels})    
    
    # fix some stuff -----
    if Experiment.get('Energy_Axis', '') == 'KineticEnergy': Experiment.update({'Energy_Axis': 'Kinetic energy (eV)'})
    if Experiment.get('Count_Rate', '') == 'Counts/Second': Experiment.update({'Count_Rate': 'Intensity (counts/s)'})
    # ----
    
    DC.update({'File': File})
    DC.update({'Experiment': Experiment})
    DC.update({'Parameters': Parameters})

    ParVals = []
    NonEnergyOrdinate = []
    for i in range(len(Parameters)): ParVals.append([])
    DATA = []
    Data = []
    with open(File.get('file'), 'rb') as f:
        #if nshup: print('reading data...')
        for i, row_ in enumerate(f):
            row = row_.decode(errors='ignore')
            row = row.replace('\r', '').replace('\n', '')

            if i>= DataRowStart:

                if row.startswith('# Cycle') and not 'Curve' in row:
                    Data = np.transpose(Data)
                    DATA.append(Data)
                    Data = []

                if row.startswith("# Parameter") and not 'Step' in row:
                    par = row.split(":")[1].split('=')[0].replace('" ', '').replace(' "', '')
                    val = row.split("=")[1]
                    if type(val) is str:
                        try: val = float(val)
                        except: pass
                    if par in Parameters:
                        ParVals[Parameters.index(par)].append(val)
                
                if row.startswith("# NonEnergyOrdinate"):
                    val = float(row.split(":")[1])
                    NonEnergyOrdinate.append(val)
                
                if not row.startswith('#'):
                    values = row.split('  ')
                    try:
                        #Data.append([float(values[0]), float(values[-1])])
                        Data.append([float(values[0]), float(values[1])])
                    except:
                        #print("There might be an issue with the number of columns in the file. Not with the")
                        #print("data file itself but how data is stored.")
                        pass
        DATA.append(np.transpose(Data))
    DATA.pop(0)
    
    DC.update({'Parameter_values': ParVals})
    NonEnergyOrdinate = np.array(NonEnergyOrdinate)
    NonEnergyOrdinate = np.unique(NonEnergyOrdinate)
    if len(NonEnergyOrdinate) > 0:
        DC.update({'Non_Energy_Ordinate': NonEnergyOrdinate})
    DC.update({'Data': DATA})

    # ------ parameters and data are now loaded. try to identify common measurements -----------------------------

    PARAMETERS = DC.get('Parameters')
    PARAMETERVALUES = DC.get('Parameter_values')
    ANALYZER = DC.get('Experiment').get('Analyzer')

    elabel = DC.get('Experiment').get('Energy_Axis')
    Meta = {}

    Measurement_type = MeasurementType(DC, shup = True)

    global defax_energy , defax_angle, defax_xdefl, defax_ydefl, defax_time, defax_position, defax_scattering_energy

    # ARPES or XPS
    sum_xps = kwargs.get("sum_xps", False)
    if Measurement_type in ['XPS', 'ARPES']:
        # ARPES: [y:e]
        data_ = []
        angle_y = np.array(DC.get('Non_Energy_Ordinate'))
        energy = np.unique(DC.get('Data')[0][0])
        NA = len(angle_y)
        NE = len(energy)
        for na in range(NA):
            data_.append(np.array(DC.get('Data')[0][1][na*NE:(na+1)*NE]))
        data = np.array(data_)
        if abs(angle_y.min()) < 1 and abs(angle_y.max()) < 1 and sum_xps: 
            DC.update({'Measurement_type': Measurement_type})
            DC.update({'x': energy})
            Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
            Meta.update({'x_type': defax_energy})
            DC.update({'int': data[:,:].sum(axis = 0)})
            Meta.update({'intlabel': DC.get('Experiment').get('Count_Rate', '')})
        elif abs(angle_y.min()) < 1 and abs(angle_y.max()) < 1 and not sum_xps:
            DC.update({'Measurement_type': Measurement_type})
            DC.update({'x': energy})
            Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
            Meta.update({'x_type': defax_energy})
            DC.update({'y': NonEnergyOrdinate})
            Meta.update({'y_label': 'Position (mm)'})
            Meta.update({'y_type': defax_position})
            DC.update({'int': data})
            Meta.update({'int_label': DC.get('Experiment').get('Count_Rate', '')})
        else:
            DC.update({'Measurement_type': Measurement_type})
            DC.update({'x': energy})
            Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
            Meta.update({'x_type': defax_energy})
            DC.update({'y': angle_y})
            Meta.update({'y_label': 'Angle (deg.)'})
            Meta.update({'y_type': defax_angle})
            DC.update({'int': data})
            Meta.update({'int_label': DC.get('Experiment').get('Count_Rate', '')})
        DC.update({'Meta': Meta})
    

    # Scattering intensity scan 
    if Measurement_type in ['scattering_intensity_scan']:
        DC.update({'Measurement_type': Measurement_type})
        timE = DC.get('Data')[0][0] * DC.get('Experiment', {}).get('Dwell_Time', 0)
        intensity = DC.get('Data')[0][1]
        DC.update({'x': timE})
        DC.update({'int': intensity})
        Meta.update({'x_label': 'Time (s)'})
        Meta.update({'x_type': defax_time})
        Meta.update({'intlabel': DC.get('Experiment').get('Count_Rate', '')})
        DC.update({'Meta': Meta})
    
    
    # Target scattering spectrum
    if Measurement_type in ['target_scattering_spectrum']:
        DC.update({'Measurement_type': Measurement_type})
        intensity = []
        for i, r in enumerate(DC.get('Data')):
            intensity.append(r[1].sum())
        intensity = np.array(intensity)
        indx = np.where(np.array(PARAMETERS) == 'ScatteringEnergy [V]')[0][0]
        scattering_energy = np.array(PARAMETERVALUES[indx])
        if scattering_energy[0] > scattering_energy[1]:
            scattering_energy = np.flip(scattering_energy)
            intensity = np.flip(intensity)
        DC.update({'x': scattering_energy})
        Meta.update({'x_label': 'Scattering energy (eV)'})
        Meta.update({'x_type': defax_scattering_energy})
        DC.update({'int': intensity})
        Meta.update({'intlabel': DC.get('Experiment').get('Count_Rate', '')})
        DC.update({'Meta': Meta})
    
    
    # Deflector scattering scan / Deflector Phoibos intensity scan 
    if Measurement_type in ['deflector_scattering_intensity_scan', 'deflector_total_intensity_scan']:
        DC.update({'Measurement_type': Measurement_type})
        Data = DC.get('Data')
        energy = np.array(np.unique(Data[0][0]))
        deflector = np.array(PARAMETERVALUES[0])
        if 'X' in PARAMETERS[0]:
            deflector_label = 'X-Deflection (deg.)'
            Meta.update({'y_type': defax_xdefl})
        else: 
            deflector_label = 'Y-Deflection (deg.)'
            Meta.update({'y_type': defax_ydefl})
        if 'ShiftY' in PARAMETERS[0]:
            deflector = -deflector                           ### SALY and ShiftY ###
        Map = []
        for d in range(len(Data)):
            Map.append(np.array(Data[d][1]))
        #
        if deflector[0] > deflector[1]:
            deflector = np.array(list(reversed(deflector)))
            Map = list(reversed(Map))
        #
        Map = np.array(Map)
        DC.update({'x': energy})
        Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
        DC.update({'y': deflector})
        Meta.update({'y_label': deflector_label})
        DC.update({'int': Map})
        Meta.update({'int_label': DC.get('Experiment').get('Count_Rate', '')})
        DC.update({'Meta': Meta})
    
    # Fermi map
    if Measurement_type == 'fermi_map':
        DC.update({'Measurement_type': Measurement_type})
        Data = DC.get('Data')
        energy = np.array(np.unique(Data[0][0]))
        Lene = len(energy)
        yAngle = DC.get('Non_Energy_Ordinate')
        Lya = len(yAngle)
        xAngle = np.array(DC.get('Parameter_values')[0])
        MAP = []
        for iShift, Shift in enumerate(xAngle):
            Map = []
            for iY in range(Lya):
                Map.append(np.array(DC.get('Data')[iShift][1][iY*Lene:(iY+1)*Lene]))
            MAP.append( np.array(Map) )
        #
        if xAngle[0] > xAngle[1]:
            xAngle = np.array(list(reversed(xAngle)))
            MAP = list(reversed(MAP))
        #
        MAP = np.array(MAP)
        DC.update({'x': energy})
        Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
        Meta.update({'x_type': defax_energy})
        DC.update({'y': yAngle})
        Meta.update({'y_label': 'Angle (deg.)'})
        Meta.update({'y_type': defax_angle})
        DC.update({'z': xAngle})
        Meta.update({'z_label': 'X-Deflection (deg.)'})
        Meta.update({'z_type': defax_xdefl})
        DC.update({'int': MAP})
        Meta.update({'int_label': DC.get('Experiment').get('Count_Rate', '')})
        DC.update({'Meta': Meta})
    
    # -------- manipulator stuff

    if Measurement_type in ['manipulator_line_scan']:         #*REWRITTEN*
        # have to extend this later to accommodate more options
        if not DC.get('Experiment', {}).get('Scan_Mode', '') == 'SnapshotFAT':
            print('Have not finished this method yet. Only works for SnapshotFAT')
            print('under certain contidions (non energy channels = 1,...)')
            print('Returning raw data. Sort it yourself.')
            return DC
        #
        DC.update({'Measurement_type': Measurement_type})
        for i, position in enumerate(PARAMETERVALUES[0]):
            energy = DC.get('Data')[0][0]
            position = np.array(DC.get('Parameter_values')[0])
            imap = []
            for d in DC.get('Data'): imap.append(d[1])
            imap = np.array(imap)
        DC.update({'Measurement_type': Measurement_type})
        DC.update({'x': energy})
        Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
        Meta.update({'x_type': defax_energy})
        DC.update({'y': position})
        Meta.update({'y_label': 'Position (mm)'})
        Meta.update({'y_type': defax_position})
        DC.update({'int': imap})
        Meta.update({'int_label': DC.get('Experiment').get('Count_Rate', '')})
        DC.update({'Meta': Meta})

    # -------- spin stuff below

    # Spin EDC (simple)  
    if Measurement_type in ['spin_edc']:
        DC.update({'Measurement_type': Measurement_type})
        energy = np.array(DC.get('Data')[0][0])
        edc_off = []; edc_on = []
        for i, polarity in enumerate(PARAMETERVALUES[0]):
            edc = DC.get('Data')[i][1]
            if 'OFF' in polarity:
                edc_off.append(edc)
            else:
                edc_on.append(edc)
        edc_off = np.array(edc_off)
        edc_on = np.array(edc_on)
        
        edc_off_avg = np.copy(edc_off[0]) * 0
        for edc in edc_off:
            edc_off_avg += edc
        edc_off_avg = edc_off_avg / len(edc_off)
        edc_on_avg = np.copy(edc_on[0]) * 0
        for edc in edc_on:
            edc_on_avg += edc
        edc_on_avg = edc_on_avg / len(edc_on)
        #
        edc_off_avg = np.array(edc_off_avg)
        edc_on_avg = np.array(edc_on_avg)
        # 
        doFilter = kwargs.get('remove_spikes', None)
        if type(doFilter) is type(None):
            doFilter = kwargs.get('filter_outliers', None)
            if not type(doFilter) is type(None):
                print(Fore.LIGHTRED_EX + "(Note: in the future use arg. remove_spikes instead of filter_outliers.)"+ Fore.RESET)
        if doFilter:
            threshold = kwargs.get('threshold', 1.2)
            edc_off_avg = removeSpikes(edc_off_avg, threshold)
            edc_on_avg = removeSpikes(edc_on_avg, threshold)
        #
        asym = (edc_off_avg - edc_on_avg) / (edc_off_avg + edc_on_avg)
        #
        DC.update({'x': energy})
        DC.update({'int': edc_off_avg})
        DC.update({'int_on': edc_on_avg})
        DC.update({'int_all': edc_off})
        DC.update({'int_all_on': edc_on})
        DC.update({'asymmetry': asym})
        Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
        Meta.update({'x_type': defax_energy})
        Meta.update({'int_label': DC.get('Experiment').get('Count_Rate', '')})
        DC.update({'Meta': Meta})


    # spin_deflector_1a
    if Measurement_type in ['spin_deflector_1a']:
        DC.update({'Measurement_type': Measurement_type})
        indxP = DC.get('Parameters').index('NegativePolarity')
        PolarityValues = DC.get('Parameter_values')[indxP]

        # 'SAL X [deg]', 'ShiftX [a.u.]', 'SAL Y [deg]',  or 'ShiftY [a.u.]
        indxD = -1
        Deflector = ''
        try:
            indxD = DC.get('Parameters').index('SAL X [deg]')
        except:
            try:
                indxD = DC.get('Parameters').index('SAL Y [deg]')
            except:
                try:
                    indxD = DC.get('Parameters').index('ShiftX [a.u.]')
                except:
                    try:
                        indxD = DC.get('Parameters').index('ShiftY [a.u.]')
                    except:
                        pass
        DeflectorValues = DC.get('Parameter_values')[indxD]

        if 'ShiftY' in DC.get('Parameters'): DeflectorValues = -DeflectorValues
        DeflectorValuesU = np.unique(DeflectorValues)
        
        energy = np.array(DC.get('Data')[0][0])
        if len(energy) == 1: energy = np.array([DC.get('Experiment', {}).get('Ek', 0)])
        edc_off = []
        edc_on = []
        for i in range(len(DeflectorValuesU)):
            edc_off.append([]); edc_on.append([])
        for i, row in enumerate(DC.get('Data')):
            deflector_index = abs(DeflectorValuesU - DeflectorValues[i]).argmin()
            if 'OFF' in PolarityValues[i]:
                edc_off[deflector_index].append(np.array(row[1]))
            else:
                edc_on[deflector_index].append(np.array(row[1]))
        
        edc_off_avg, edc_on_avg = [], []
        for i in range(len(DeflectorValuesU)):
            edc_avg = np.copy(energy) * 0
            cnt = 0
            for edc in edc_off[i]:
                edc_avg += edc
                cnt += 1
            edc_off_avg.append(edc_avg/cnt)
            edc_avg = np.copy(energy) * 0
            cnt = 0
            for edc in edc_on[i]:
                edc_avg += edc
                cnt += 1
            edc_on_avg.append(edc_avg/cnt)
        edc_off_avg = np.array(edc_off_avg)
        edc_on_avg = np.array(edc_on_avg)
        #
        if len(energy) == 1:
            edc_off_avg = edc_off_avg.flatten()
            edc_on_avg = edc_on_avg.flatten()
        #
        doFilter = kwargs.get('remove_spikes', None)
        if type(doFilter) is type(None):
            doFilter = kwargs.get('filter_outliers', None)
            if not type(doFilter) is type(None):
                print(Fore.LIGHTRED_EX + "(Note: in the future use arg. remove_spikes instead of filter_outliers.)"+ Fore.RESET)
        if doFilter:
            threshold = kwargs.get('threshold', 1.2)
            edc_off_avg = filterOutliers(edc_off_avg, threshold)
            edc_on_avg = filterOutliers(edc_on_avg, threshold)
        #
        asym = (edc_off_avg - edc_on_avg) / (edc_off_avg + edc_on_avg)
        #
        DC.update({'x': np.array(energy)})
        Meta.update({'x_type': defax_energy})
        if 'X' in DC.get('Parameters')[indxD]:
            ShiftParameterLabel = 'X-Deflection (deg.)'
            Meta.update({'y_type': defax_xdefl})
        else: 
            ShiftParameterLabel = 'Y-Deflection (deg.)'
            Meta.update({'y_type': defax_ydefl})
        DC.update({'y': DeflectorValuesU})
        DC.update({'int': edc_off_avg})
        DC.update({'int_on': edc_on_avg})
        DC.update({'int_all': edc_off})
        DC.update({'int_all_on': edc_on})
        DC.update({'asymmetry': asym})
        Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
        Meta.update({'y_label': ShiftParameterLabel})
        Meta.update({'int_label': DC.get('Experiment').get('Count_Rate', '')})
        DC.update({'Meta': Meta})
    
    # spin_deflector_1b
    if Measurement_type in ['spin_deflector_1b']:
        #print("\nNOTE THAT loading of 2d deflector maps are is implemented yet, only deflector line scans (w. 2 deflectors).\n")

        DC.update({'Measurement_type': Measurement_type})
        indx = DC.get('Parameters').index('NegativePolarity')
        PolarityValues_ = DC.get('Parameter_values')[indx]
        indx = DC.get('Parameters').index('ShiftX [a.u.]')
        ShiftX_ = DC.get('Parameter_values')[indx]
        indx = DC.get('Parameters').index('ShiftY [a.u.]')
        ShiftY_ = DC.get('Parameter_values')[indx]
        #
        PolarityValues = np.unique(PolarityValues_)
        ShiftX = np.array(np.unique(ShiftX_))
        ShiftY = np.array(np.unique(ShiftY_))
        if ShiftX_[0] > ShiftX_[-1]: ShiftX = ShiftX[::-1]
        if ShiftY_[0] > ShiftY_[-1]: ShiftY = ShiftY[::-1]

        energy = np.array(DC.get('Data')[0][0])
        if len(energy) == 1: energy = [np.array(DC.get('Experiment', {}).get('Ek', 0))]
        edc_off = []
        edc_on = []
        for i in range(len(ShiftX)):
            edc_off.append([]); edc_on.append([])
        for i, row in enumerate(DC.get('Data')):
            deflector_index = abs(ShiftX - ShiftX_[i]).argmin()
            if 'OFF' in PolarityValues_[i]:
                edc_off[deflector_index].append(np.array(row[1]))
            else:
                edc_on[deflector_index].append(np.array(row[1]))
        #
        edc_off_avg, edc_on_avg = [], []
        for i in range(len(ShiftX)):
            edc_avg = np.copy(energy) * 0
            cnt = 0
            for edc in edc_off[i]:
                edc_avg += edc
                cnt += 1
            edc_off_avg.append(edc_avg/cnt)
            edc_avg = np.copy(energy) * 0
            cnt = 0
            for edc in edc_on[i]:
                edc_avg += edc
                cnt += 1
            edc_on_avg.append(edc_avg/cnt)
        edc_off_avg = np.array(edc_off_avg)
        edc_on_avg = np.array(edc_on_avg)
        #
        if len(energy) == 1:
            edc_off_avg = edc_off_avg.flatten()
            edc_on_avg = edc_on_avg.flatten()
        #
        doFilter = kwargs.get('remove_spikes', None)
        if type(doFilter) is type(None):
            doFilter = kwargs.get('filter_outliers', None) # legacy argument
            if not type(doFilter) is type(None):
                print(Fore.LIGHTRED_EX + "(Note: in the future use arg. remove_spikes instead of filter_outliers.)"+ Fore.RESET)

        if doFilter:
            threshold = kwargs.get('threshold', 1.2)
            edc_off_avg = removeSpikes(edc_off_avg, threshold)
            edc_on_avg = removeSpikes(edc_on_avg, threshold)
        #
        asym = (edc_off_avg - edc_on_avg) / (edc_off_avg + edc_on_avg)
        #
        #
        DC.update({'x': np.array(energy)})
        DC.update({'y': ShiftX})
        DC.update({'y2': ShiftY})
        DC.update({'int': edc_off_avg})
        DC.update({'int_on': edc_on_avg})
        DC.update({'int_all': edc_off})
        DC.update({'int_all_on': edc_on})
        DC.update({'asymmetry': asym})
        Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
        Meta.update({'x_type': defax_energy})
        Meta.update({'y_label': 'X-Deflection (deg.)'})
        Meta.update({'y_type': defax_xdefl})
        Meta.update({'y2_label': 'Y-Deflection (deg.)'})
        Meta.update({'y2_type': defax_xdefl})
        Meta.update({'int_label': DC.get('Experiment').get('Count_Rate', '')})
        DC.update({'Meta': Meta})
    
    DC.update({'Type': DC.get('Measurement_type')})

    if nshup:
        Info(DC)
    
    return DC



# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================

def info(D = {}):
    """"""
    RED = Fore.RED; BLU = Fore.BLUE; BLK = Fore.RESET
    if not type(D) is dict: 
        print(f'{RED}Argument D must be a grumpy dict.{BLK}'); return
    if D.get('Experiment', {}) == {}:
        print(f"{RED}Argument D is probably not a grumpy dict but let's check it anyway...{BLK}")
    print(f"{BLK}Type:  {BLU}{D.get('Type', '?')}{BLK}")
    ax = _getAxes(D = D)
    print('Axes:')
    if len(ax) == 0:
        print(f'{RED}  none{BLK}')
    else:
        for a in ax:
            print("{0}  {1}: {2} points, from {3} to {4}, {5}{6}".format(BLU, a[0], a[1], a[2], a[3], a[4], BLK))
    print('Intensity:')
    I = D.get('int')
    if type(I) is type(None):
        print(f'{RED}  none{BLK}')
    else:
        print(f"{BLU}  int: {np.shape(I)}, min = {np.nanmin(I)}, max = {np.nanmax(I)}, {D.get('Meta',{}).get('int_label', '')}{BLK}")
    if D.get('Type', '?').lower().startswith('spin'):
        print('Spin intensities:')
        for key in ['int', 'int_on', 'int_all', 'int_all_on', 'asymmetry']:
            print(BLU + "{0}  {1:<10}: {2}{3}".format(BLU, key, np.shape(D.get(key,[]), BLK)))



def Info(D = {}):
    """"""
    RED = Fore.RED
    BLK = Fore.RESET
    BLU = Fore.BLUE
    if not type(D) is dict: 
        print(RED + 'Argument D must be a grumpy dict.'+ BLK); return
    #
    Experiment = D.get('Experiment', {})
    if Experiment == {}:
        print(RED + 'This is probably not a grumpy dict.' + BLK); return
    File = D.get('File', {})
    #
    print('----------------')
    print('Prodigy version:    {0}{1}{2}'.format(BLU, Experiment.get('Version', '?'), BLK))
    print('File:               {0}{1}{2}'.format(BLU, File.get('file_name', '?'), BLK))
    print('Spectrum ID:        {0}{1}{2}'.format(BLU, Experiment.get('Spectrum_ID', '?'), BLK))
    print('Measurement type:   {0}The data was identified as {1}{2}'.format(BLU, D.get('Measurement_type', '?'), BLK))
    if not D.get('Measurement_type', '?') == D.get('Type', '??'):
        print('Type:               {0}The data is now of type {1}{2}'.format(BLU, D.get('Type', '?'), BLK)) 
    print()
    print('Analyzer:           {0}{1}{2}'.format(BLU, Experiment.get('Analyzer', '?'), BLK))
    print('Lens mode:          {0}{1}{2}'.format(BLU, Experiment.get('Lens_Mode', '?'), BLK))
    print('Scan mode:          {0}{1}{2}'.format(BLU, Experiment.get('Scan_Mode', '?'), BLK))
    print('Energy              {0}{1}{2}'.format(BLU, Experiment.get('Energy_Axis', '?').strip('(eV)'), BLK))
    print('Photon energy:      {0}{1} eV{2}'.format(BLU, Experiment.get('Excitation_Energy', '?'), BLK))
    print('Kinetic energy:     {0}{1} eV{2}'.format(BLU, Experiment.get('Kinetic_Energy', '?'), BLK))
    print('Pass energy:        {0}{1} eV{2}'.format(BLU, Experiment.get('Pass_Energy', '?'), BLK))
    print()
    print('Num. of parameters: {0}{1}{2}'.format(BLU, len(D.get('Parameters', [])), BLK))
    if len(D.get('Parameters', [])) > 0:
        print('Parameters:         {0}{1}{2}'.format(BLU, ', '.join(p for p in D.get('Parameters', []) ), BLK))
    print('\nData:' + BLU)
    print('Loaded data:        {0}'.format(np.shape(D.get('Data', []))))
    if not np.shape(D.get('x')) == ():
        print('axis x:             {0}, {1}, {2} - {3}, {4}'.format(np.shape(D.get('x', [])), D.get('Meta', {}).get('x_type', '?'), D.get('x').min(), D.get('x').max(), D.get('Meta', {}).get('x_label', '?')))
    if not np.shape(D.get('y')) == ():
        print('axis y:             {0}, {1}, {2} - {3}, {4}'.format(np.shape(D.get('y', [])), D.get('Meta', {}).get('y_type', '?'), D.get('y').min(), D.get('y').max(), D.get('Meta', {}).get('y_label', '?')))
    if not np.shape(D.get('y2')) == ():
        print('axis y2:            {0}, {1}, {2} - {3}, {4}'.format(np.shape(D.get('y2', [])), D.get('Meta', {}).get('y2_type', '?'), D.get('y2').min(), D.get('y2').max(), D.get('Meta', {}).get('y2_label', '?')))
    if not np.shape(D.get('z')) == ():
        print('axis z:             {0}, {1}, {2} - {3}, {4}'.format(np.shape(D.get('z', [])), D.get('Meta', {}).get('z_type', '?'), D.get('z').min(), D.get('z').max(), D.get('Meta', {}).get('z_label', '?')))
    if not np.shape(D.get('int')) == ():
        shp = np.shape(D.get('int'))
    else:
        shp = "value = {0}".format(D.get('int'))
    print('int:                {0}'.format(shp) + Fore.RESET)
    print(BLK)

    

# ============================================================================================================


def _getAxes(D = {}):
    """
    Helper
    """
    axes = []
    M = D.get('Meta', {})
    if 'x' in D: axes.append(['x', len(D['x']), D['x'].min(), D['x'].max(), M.get('x_label', '')])
    if 'y' in D: axes.append(['y', len(D['y']), D['y'].min(), D['y'].max(), M.get('y_label', '')])
    if 'z' in D: axes.append(['z', len(D['z']), D['z'].min(), D['z'].max(), M.get('z_label', '')])
    return axes


def _fixIndices(axis, x1, x2, ix1, ix2):
    """
    Helper.
    """
    NoneType = type(None)
    if not type(x1) is NoneType or not type(x2) is NoneType:
        if type(x1) is NoneType: x1 = axis.min()
        if type(x2) is NoneType: x2 = axis.max()
        if x1 > x2: x1, x2 = x2, x1
        if x1 < axis.min(): x1 = axis.min()
        if x2 > axis.max(): x2 = axis.max()
        ix1 = abs(axis - x1).argmin()
        ix2 = abs(axis - x2).argmin()
        if ix1 > ix2: ix1, ix2 = ix2, ix1
    elif not type(ix1) is NoneType or not type(ix2) is NoneType:
        if type(ix1) is NoneType: ix1 = 0
        if type(ix2) is NoneType: ix2 = len(axis) - 1
        if ix1 > ix2: ix1, ix2 = ix2, ix1
        x1 = axis[ix1]; x2 = axis[ix2]
        if x1 > x2: x1, x2 = x2, x1
    else:
        x1 = axis.min(); x2 = axis.max()
        ix1 = 0; ix2 = len(axis) - 1
    return x1, x2, ix1, ix2



# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================


def SubArray(D = {}, shup = False, **kwargs):
    """
    Return a subset of an array, defined by x1 and x2 (1d,2d,3d), and y1 and y2 (2d,3d), and z1 and z2 (3d).
    Instead of values x1, x2, y1, ..., z2 you can use indices ix1, ix2, iy1, ..., iz2.
    """
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must be a grumpy dict.' + Fore.RESET); return D
    #
    av_axes = _getAxes(D)
    if len(av_axes) == 0:
        print(Fore.RED + 'Could not find any axes in the dict.' + Fore.RESET); return D
    #
    recognized_kwargs = ['x1', 'x2', 'ix1', 'ix2', 'y1', 'y2', 'iy1', 'iy2', 'z1', 'z2', 'iz1', 'iz2']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    x1 = kwargs.get('x1');      x2 = kwargs.get('x2')
    ix1 = kwargs.get('ix1');    ix2 = kwargs.get('ix2')
    y1 = kwargs.get('y1');      y2 = kwargs.get('y2')
    iy1 = kwargs.get('iy1');    iy2 = kwargs.get('iy2')
    z1 = kwargs.get('z1');      z2 = kwargs.get('z2')
    iz1 = kwargs.get('iz1');    iz2 = kwargs.get('iz2')
    #
    if type(D.get('int')) is None:
        print(Fore.RED + "There is no intensity data in the dict." + Fore.RESET); return D
    # 
    DD = copy.deepcopy(D)
    #
    #_fixIndices
    if len(np.shape(DD['int'])) > 0:
        x1, x2, ix1, ix2 = _fixIndices(DD['x'], x1, x2, ix1, ix2)
        DD.update({'x': DD['x'][ix1:ix2+1]})
        if not shup: print('x-range: {0} to {1}'.format(x1, x2))
    if len(np.shape(DD['int'])) > 1:
        y1, y2, iy1, iy2 = _fixIndices(DD['y'], y1, y2, iy1, iy2)
        DD.update({'y': DD['y'][iy1:iy2+1]})
        if not shup: print('y-range: {0} to {1}'.format(y1, y2))
    if len(np.shape(DD['int'])) > 2:
        z1, z2, iz1, iz2 = _fixIndices(DD['z'], z1, z2, iz1, iz2)
        DD.update({'z': DD['z'][iz1:iz2+1]})
        if not shup: print('z-range: {0} to {1}'.format(DD['z'].min(), DD['z'].max()))
    #
    # Non-spin data
    if not DD.get('Type', '').lower().startswith('spin'):
        if len(np.shape(DD['int'])) == 1:
            DD.update({'int': DD['int'][ix1:ix2+1]}) 
        elif len(np.shape(DD['int'])) == 2:
            DD.update({'int': DD['int'][iy1:iy2+1, ix1:ix2+1]})
        elif len(np.shape(DD['int'])) == 3:
            DD.update({'int': DD['int'][iz1:iz2+1, iy1:iy2+1, ix1:ix2+1]})
    # spin data
    else:
        # simple spin edc
        if len(np.shape(DD['int'])) == 1:
            DD.update({'int': DD.get('int')[ix1:ix2+1]})
            DD.update({'int_on': DD.get('int_on')[ix1:ix2+1]})
            DD.update({'int_all': DD.get('int_all')[:,ix1:ix2+1]})
            DD.update({'int_all_on': DD.get('int_all_on')[:,ix1:ix2+1]})
            DD.update({'asymmetry': DD.get('asymmetry')[ix1:ix2+1]})
        # deflector spin
        else:
            print(Fore.RED + 'Deflector spin not implemented yet.' + Fore.RESET); return D
            
    #
    if not shup:
        print('int:', np.shape(DD['int']))
    return DD



# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================



def Compact(D = {}, axis = '', shup = False, **kwargs):   # rewritten
    """
    Compact an array along an axis defined by argument axis. For 1d data the argument is x (used by default),
    for 2d data the axis can be x or y, and for 3d data axis can be x, y, or z.
    """
    recognized_kwargs = ['']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    if not type(D) is dict:
        print(Fore.RED + 'Argument D must be a grumpy dict.'+ Fore.RESET); return D
    #
    av_axes = _getAxes(D)
    if len(av_axes) == 0:
        print(Fore.RED + 'Could not find any axes in the dict.' + Fore.RESET); return D
    #
    if D.get('Type').startswith('spin'):
        print(Fore.RED + 'Compact does not work with spin data yet.' + Fore.RESET); return D
    #
    axis = axis.lower()
    if not axis in list(np.transpose(av_axes)[0]):
        print(Fore.RED + "Axis {0} is not found in the dict.".format(axis) + Fore.RESET); return D
    #
    Int = D.get('int')
    if type(Int) is None:
        print(Fore.RED + "There is no intensity data in the dict." + Fore.RESET); return D
    # 
    DD = copy.deepcopy(D)
    Meta = copy.deepcopy(DD.get('Meta', {}))
    #
    if len(np.shape(D.get('int'))) == 1:
        DD.update({'int': D.get('int').sum(axis = 0)})
        del DD['x']
        del Meta['x_label']
        del Meta['x_type']
        DD.update({'Meta': Meta})
    #
    elif len(np.shape(D.get('int'))) == 2:
        if axis == 'x':
            DD.update({'x': D.get('y')})
            Meta.update({'x_label': D.get('Meta', {}).get('y_label', '')})
            Meta.update({'x_type': D.get('Meta', {}).get('y_type', '')})
            DD.update({'int': D.get('int').sum(axis = 1)})
            del DD['y']
            del Meta['y_label']
            del Meta['y_type']
            DD.update({'Meta': Meta})

        else: #(axis == y)
            DD.update({'x': D.get('x')})
            Meta.update({'x_label': D.get('Meta', {}).get('x_label', '')})
            Meta.update({'x_type': D.get('Meta', {}).get('x_type', '')})
            DD.update({'int': D.get('int').sum(axis = 0)})
            del DD['y']
            del Meta['y_label']
            del Meta['y_type']
            DD.update({'Meta': Meta})
    #
    elif len(np.shape(D.get('int'))) == 3:      # CORRECT TYPE UPDATING
        if axis == 'x':
            DD.update({'x': DD.get('y')})
            DD.update({'y': DD.get('z')})
            Meta.update({'x_label': DD.get('Meta', {}).get('y_label', '')})
            Meta.update({'x_type': DD.get('Meta', {}).get('y_type', '')})
            Meta.update({'y_label': DD.get('Meta', {}).get('z_label', '')})
            Meta.update({'y_type': DD.get('Meta', {}).get('z_type', '')})
            DD.update({'int': DD.get('int').sum(axis = 2)})
        elif axis == 'y':
            DD.update({'x': DD.get('x')})
            Meta.update({'x_label': DD.get('Meta', {}).get('x_label', '')})
            Meta.update({'x_type': DD.get('Meta', {}).get('x_type', '')})
            DD.update({'y': DD.get('z')})
            Meta.update({'y_label': DD.get('Meta', {}).get('z_label', '')})
            Meta.update({'y_type': DD.get('Meta', {}).get('z_type', '')})
            DD.update({'int': DD.get('int').sum(axis = 1)})
        elif axis == 'z':
            DD.update({'x': DD.get('x')})
            Meta.update({'x_label': DD.get('Meta', {}).get('x_label', '')})
            Meta.update({'x_type': DD.get('Meta', {}).get('x_type', '')})
            DD.update({'y': DD.get('y')})
            Meta.update({'y_label': DD.get('Meta', {}).get('y_label', '')})
            Meta.update({'y_type': DD.get('Meta', {}).get('y_type', '')})
            DD.update({'int': DD.get('int').sum(axis = 0)})
        del DD['z']
        del Meta['z_label']
        del Meta['z_type']
        DD.update({'Meta': Meta})
    #
    Type = DD.get('Type')
    if Type in ['XPS', 'scattering_intensity_scan', 'target_scattering_spectrum']: 
        DD.update({'Type': 'value'})
    elif Type in ['ARPES', 'fermi_map_slice']: 
        DD.update({'Type': 'profile'})
    elif Type == 'fermi_map':
        DD.update({'Type': 'fermi_map_slice'})

    #
    return DD






# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================


def Profile(D = {}, **kwargs):
    """
    arguments:
    D       grumpy dict
    axis    compacting this axis
    x1, x2  x-limits, optional
    y1, y2  y-limits, optional
    z1, z2  z-limits, optional
    """
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must be a grumpy dict.' + Fore.RESET); return D
    #
    recognized_kwargs = ['axis', 'x1', 'x2', 'y1', 'y2', 'z1', 'z2']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    axis = kwargs.get('axis', '')
    x1, x2 = kwargs.get('x1'), kwargs.get('x2')
    y1, y2 = kwargs.get('y1'), kwargs.get('y2')
    z1, z2 = kwargs.get('z1'), kwargs.get('z2')
    #
    axes_info = _getAxes(D = D)
    if len(axes_info) == 0:
        print(Fore.RED + 'Can not find any axes in argument D.' + Fore.RESET); return D
    axes = list(np.transpose(axes_info))[0]
    #
    if axis == '':
        print(Fore.RED + 'Argument axis must be one of the following: {0}'.format(axes) + Fore.RESET); return D
    #
    if not axis in axes:
        print(Fore.RED + 'Argument axis must be one of the following: {0}'.format(axes) + Fore.RESET); return D
    #
    ntp = type(None)
    DD = copy.deepcopy(D)
    if len(axes) > 2:
        if not type(z1) is ntp or not type(z2) is ntp:
            DD = SubArray(D = DD, z1 = z1, z2 = z2, shup = True)
    if len(axes) > 1:
        if not type(y1) is ntp or not type(y2) is ntp:
            DD = SubArray(D = DD, y1 = y1, y2 = y2, shup = True)
    if len(axes) > 0:
        if not type(x1) is ntp or not type(x2) is ntp:
            DD = SubArray(D = DD, x1 = x1, x2 = x2, shup = True)
    #
    DD = Compact(D = DD, axis = axis)    

    return DD






# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================


def ShiftAxis(D = {}, axis = '', v = 0, shup = False, **kwargs):
    """
    Shifts the values on the axis with v.
    """
    recognized_kwargs = ['']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    if not type(D) is dict:
        print(Fore.RED + 'Argument D must be a grumpy dict.' + Fore.RESET); return D
    #
    axes_info = _getAxes(D = D)
    axes = list(np.transpose(axes_info))[0]
    if not axis in axes:
        print(Fore.RED + "axis '{0}' is not found in this grumpy dict.".format(axis) + Fore.RESET); return D
    #
    B = copy.deepcopy(D)
    B.update({axis: np.array(B.get(axis)) + v})
    if not shup:
        print('axis {0} is shifted by {1}, from starting and ending at {2} and {3} to'.format(axis, v, D.get(axis)[0], D.get(axis)[-1]))
        print('starting and ending at {0} and {1}.'.format(B.get(axis)[0], B.get(axis)[-1]))
    return B


# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================



def Plot(D = {}, ax = None, shup = False, **kwargs):
    """
    To quickly plot data in a dict from Load() or from dicts by e.g. Profile(), Compact(),
    QuickSpin(), ...
    Returns an ax.

    General kwargs:
        figsize         always applicable, type tuple, (x,y)
        axvline         always applicable, type list, [x1, x2,...]
        axhline         always applicable, type list, [y1, y2,...]
        title           always applicable, type str
        cmap            applicable for 2d plots, type str
        colorbar        applicable for 2d plots, type bool
        vmin            applicable for 2d plots, type number
        vmax            applicable for 2d plots, type number
        rotate          rotate 2d plot 90 deg, type bool
        intensity       type str, applicable for: 
                            spindata,               values: 'asym', 'all', 'scans', 'avg'
                            deflector_spin_edc_1,   values: 'asym', 'onoff', 'on', 'off'
                            deflector_spin_edc_2,   values: 'asym', 'onoff', 'on', 'off'
    """
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must be a grumpy dict.'+ Fore.RESET); return ax

    
    # kwargs ----------
    recognized_kwargs = ['intensity', 'cmap', 'colorbar', 'vmin', 'vmax', 'rotate', 'aspect',
                         'color', 'linewidth', 'linestyle', 'axvline', 'axhline', 'title', 'figsize']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    kw_intensity = kwargs.get('intensity', '') 
    cmap = kwargs.get('cmap', 'bone_r')
    colorbar = kwargs.get('colorbar', False)
    vmin = kwargs.get('vmin', None)
    vmax = kwargs.get('vmax', None)
    rotate = kwargs.get('rotate', None)
    aspect = kwargs.get('aspect', '')
    kw_color = kwargs.get('color', 'tab:blue')
    kw_linewidth = kwargs.get('linewidth', 0.75)
    kw_linestyle = kwargs.get('linestyle', '-')

    kw_axvline = kwargs.get('axvline', [])
    kw_axhline = kwargs.get('axhline', [])
    kw_title = kwargs.get('title', '')
    # -----------------

    # if this is a fit-result dict...
    if 'yfit' in D and 'par' in D and 'method' in D:
        return PlotFitRes(D = D, ax = ax, **kwargs)

    Type = D.get('Type', '')

    if Type == 'fermi_map':
        PlotFermiMap(D = D, **kwargs)
        return None
    
    # ----------------- plot 1d data
    axes_in_D = _getAxes(D)
    if len(axes_in_D) == 1:
        if not hasattr(ax, 'plot'):
            figsize = kwargs.get('figsize', (4, 2.5))
            fig, ax = plt.subplots(figsize = figsize)

        if not Type.lower().startswith('spin'):
            ax.plot(D.get('x'), D.get('int'), linewidth = kw_linewidth, color = kw_color, linestyle = kw_linestyle)
            ax.set_xlabel(D.get('Meta',{}).get('x_label', ''))
            ax.set_ylabel(D.get('Meta',{}).get('int_label', ''))
            ax.set_title('ID: {0}'.format(D.get('Experiment').get('Spectrum_ID')))
    
        elif Type == 'spin_edc':
            if not hasattr(ax, 'plot'):
                figsize = kwargs.get('figsize', (4,2.5))
                fig, ax = plt.subplots(figsize = figsize)
            en = D.get('x')
            edc_offs = D.get('int_all')
            edc_ons = D.get('int_all_on')
            edc_off = D.get('int')
            edc_on = D.get('int_on')
            asym = D.get('asymmetry')
            if not kw_intensity in ['all', 'scans', 'avg', 'asym']: kw_intensity = 'all'
            if kw_intensity in ['all', 'scans']:
                for edc in edc_offs:
                    ax.plot(en, edc, linewidth = 0.5, color = 'tab:blue')
                for edc in edc_ons:
                    ax.plot(en, edc, linewidth = 0.5, color = 'tab:orange')
                ax.set_ylabel(D.get('Meta',{}).get('int_label', ''))
            if kw_intensity in ['all', 'avg']:
                ax.plot(en, edc_off, linewidth = 1.25, color = 'tab:blue', label = 'polarity +')
                ax.plot(en, edc_on, linewidth = 1.25, color = 'tab:orange', label = 'polarity -')
                ax.legend()
                ax.set_ylabel(D.get('Meta',{}).get('int_label', ''))
            if kw_intensity == 'asym':
                ax.plot(en, asym, linewidth = 1.25, color = kw_color)
                ax.set_ylabel('Asymmetry')
            ax.set_xlabel('Energy (eV)')
            ax.set_title('ID: {0}'.format(D.get('Experiment').get('Spectrum_ID')))

    # ----------------- plot 2d data 
    elif len(axes_in_D) == 2 and not Type.lower().startswith('spin'):
        if not hasattr(ax, 'plot'):
            figsize = kwargs.get('figsize', (4,2.5))
            fig, ax = plt.subplots(figsize = figsize)
        xaxis = D.get('x')
        yaxis = D.get('y')
        intensity = D.get('int')
        xlabel = D.get('Meta', {}).get('x_label', '')
        ylabel = D.get('Meta', {}).get('y_label', '')
        title = "ID {0}".format(D.get('Experiment', {}).get('Spectrum_ID', '-'))
        #
        global defax_angle, defax_xdefl, defax_ydefl
        rotate_ = False
        if D.get('Meta',{}).get('y_type') in [defax_angle, defax_xdefl, defax_ydefl]:
            rotate_ = True
        #
        if not type(rotate) is type(None):
            if not rotate: rotate_ = False
        if rotate_:
            xaxis, yaxis = yaxis, np.flip(xaxis)
            xlabel, ylabel = ylabel, xlabel
            intensity = np.flip(intensity.transpose(), axis = 0)
        #
        aspect_ = 'auto'
        if not aspect == '':
            aspect_ = aspect
        else:
            if D.get('Meta', {}).get('x_type', '') in [defax_angle, defax_xdefl, defax_ydefl] and \
            D.get('Meta', {}).get('y_type', '') in [defax_angle, defax_xdefl, defax_ydefl]:
                aspect_ = 'equal'
        #
        extent = [xaxis[0], xaxis[-1], yaxis[-1], yaxis[0]]
        im = ax.imshow(intensity, extent = extent, aspect = aspect_, vmin = vmin, vmax = vmax, cmap = cmap)
        if not rotate_:
            ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    

    elif (len(axes_in_D) == 2 or len(axes_in_D) == 3) and Type.lower().startswith('spin'):

        if Type in ['spin_deflector_1a', 'spin_deflector_1b']:
            if not hasattr(ax, 'plot'):
                figsize = kwargs.get('figsize', (4,2.5))
                fig, ax = plt.subplots(figsize = figsize)
            #
            energy  = D.get('x')
            defl    = D.get('y', np.array([]))
            defl_2  = D.get('y2', np.array([]))
            edc_off = D.get('int', np.array([]))
            edc_on  = D.get('int_on', np.array([]))
            asym    = D.get('asymmetry', np.array([]))
            if len(defl) > 0 and len(defl_2) == 0:
                defl_label = D.get('Meta', {}).get('y_label', '')
            else:
                defl_label = 'XY-Defleftion (deg.)'
            #
            if len(energy) == 1:   # fixed energy
                if not kw_intensity in ['asym', 'off', 'on', 'onoff', 'offon']: kw_intensity = 'onoff'
                if kw_intensity.startswith('asym'):
                    ax.plot(defl, asym, color = 'k', linewidth = kw_linewidth, label = 'Asym.')
                    ylabel = 'Asymmetry'
                if kw_intensity in ['onoff', 'offon', 'off']:
                    ax.plot(defl, edc_off, color = 'tab:blue', linewidth = kw_linewidth, label = 'Off')
                    ylabel = D.get('Meta', {}).get('int_label', '')
                if kw_intensity in ['onoff', 'offon', 'on']:
                    ax.plot(defl, edc_on, color = 'tab:red', linewidth = kw_linewidth, label = 'On')
                    ylabel = D.get('Meta', {}).get('int_label', '')
                ax.set_ylabel(ylabel)
                if defl_label.startswith('XY'):
                    if ylabel.startswith('X'):
                        xlabel = 'Deflection from ({0},{1}) to ({2},{3}) (deg.)'.format(defl[0], defl_2[0], defl[-1], defl_2[-1])
                    else:
                        xlabel = 'Deflection from ({1},{0}) to ({3},{2}) (deg.)'.format(defl[0], defl_2[0], defl[-1], defl_2[-1])
                else:
                    xlabel = D.get('Meta', {}).get('y_label', '')
                ax.set_xlabel(xlabel)
                ttl = "ID = {0}, Ek = {1}".format(D.get('Experiment',{}).get('Spectrum_ID','?'), D.get('x')[0])
                ax.set_title(ttl)
                if not ylabel == 'Asymmetry': ax.legend()
            
            else: # edc
                if not kw_intensity in ['asym', 'off', 'on']: kw_intensity = 'asym'
                intensity_map = []
                for i in range(len(asym)):
                    if kw_intensity == 'asym': curve = asym[i]
                    elif kw_intensity == 'off': curve = edc_off[i]
                    elif kw_intensity == 'on': curve = edc_off[i]
                    intensity_map.append(curve)
                intensity_map = np.array(intensity_map)
                extent = [energy[0], energy[-1], defl[-1], defl[0]]
                ims = ax.imshow(intensity_map, extent = extent, aspect = 'auto', vmin = vmin, vmax = vmax, cmap = cmap)
                if colorbar:
                    fig.colorbar(ims, ax = ax)
                ax.invert_yaxis()
                ax.set_xlabel(D.get('Meta', {}).get('x_label', ''))
                ax.set_ylabel('Deflection {0}'.format(defl_label))
                ttl = "ID = {0}".format(D.get('Experiment',{}).get('Spectrum_ID','?'))
                ax.set_title(ttl)

    elif len(axes_in_D) > 3:
        print('Too many axes for this poor method.'); return ax
        
    else:
        print(Fore.RED + 'The data is not recognized by this method. Plot it yourself.' + Fore.RESET)
    
    # do some kwarg stuff
    if not type(kw_axvline) is list:
        kw_axvline = [kw_axvline]
    if len(kw_axvline) > 0:
        for vline in kw_axvline: ax.axvline(x = vline)
    if not type(kw_axhline) is list:
        kw_axhline = [kw_axhline]
    if len(kw_axhline) > 0:
        for hline in kw_axhline: ax.axhline(y = hline)
    if not kw_title == '':
        ax.set_title(kw_title)
    
    return ax

    






# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================

def PlotFermiMap(D = {}, **kwargs):
    """
    Interactive plot of 3d data (energy, angle, x delflection).
    """
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must be of type grumpy dict.' + Fore.RESET); return
    Type = D.get('Measurement_type', '')
    if not Type == 'fermi_map':
        Type = D.get('Type', '')
        if not Type == 'fermi_map':
            print(Fore.RED + 'The data in the dict is not a Fermi map.' + Fore.RESET); return
    
    recognized_kwargs = ['figsize', 'cmap', 'linewidth']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    figsize = kwargs.get('figsize', (10,5))
    cmap = kwargs.get('cmap', 'bone_r')
    linewidth = kwargs.get('linewidth', 0.75)

    ENERGY = D.get('x')
    ANGLEX = D.get('z')
    ANGLEY = D.get('y')
    #print("\n debug -", ENERGY[0], ANGLEX[0], ANGLEY[0], "\n")     # DEBUG
    dENERGY = (ENERGY[-1]-ENERGY[0])/len(ENERGY)
    dANGLEX = abs(ANGLEX[1]-ANGLEX[0])
    dANGLEY = abs(ANGLEY[1]-ANGLEY[0])

    SliderE = ipw.FloatSlider(min=ENERGY[0], max=ENERGY[-1], step = dENERGY, description = 'Energy', value = ENERGY.mean(), readout_format = ".3f")
    SliderX = ipw.FloatSlider(min=ANGLEX.min(), max=ANGLEX.max(), step = dANGLEX, description = 'Deflection X', value = ANGLEX.mean(), readout_format = ".2f")
    SliderY = ipw.FloatSlider(min=ANGLEY[0], max=ANGLEY[-1], step = dANGLEY, description = 'Angle Y', value = ANGLEY.mean(), readout_format = ".2f")

    SliderDE = ipw.FloatSlider(min=0, max=(ENERGY[-1]-ENERGY[0]), step = dENERGY, description = 'dE', value = 5*dENERGY, readout_format = ".3f")
    SliderDX = ipw.FloatSlider(min=0, max=(ANGLEX.max()-ANGLEX.min()), step = dANGLEX, description = 'dX', value = dANGLEX, readout_format = ".2f")
    SliderDY = ipw.FloatSlider(min=0, max=(ANGLEY[-1]-ANGLEY[0]), step = dANGLEY, description = 'dX', value = 5*dANGLEY, readout_format = ".2f")

    box_sliders_x = ipw.HBox([SliderE, SliderX, SliderY])
    box_sliders_dx = ipw.HBox([SliderDE, SliderDX, SliderDY])
    box_sliders = ipw.VBox([box_sliders_x, box_sliders_dx])

    extentE = [ANGLEX[0], ANGLEX[-1], ANGLEY[-1], ANGLEY[0]] 
    extentX = [ANGLEY[0], ANGLEY[-1], ENERGY[-1], ENERGY[0]]
    extentY = [ANGLEX[0], ANGLEX[-1], ENERGY[-1], ENERGY[0]]

    def plot(X, Y, E, DX, DY, DE):
        fig, ax = plt.subplots(figsize = figsize, ncols = 3)
        plt.tight_layout()
        #
        XY = Compact(D = SubArray(D = D, x1 = E-DE/2, x2 = E+DE/2, shup = True), axis = 'x', shup = True)
        YE = Compact(D = SubArray(D = D, z1 = X-DX/2, z2 = X+DX/2, shup = True), axis = 'z', shup = True)
        XE = Compact(D = SubArray(D = D, y1 = Y-DY/2, y2 = Y+DY/2, shup = True), axis = 'y', shup = True)
        #
        _ = ax[0].imshow(XY['int'].transpose(), extent = extentE, aspect = 'equal', cmap = cmap)
        _ = ax[1].imshow(YE['int'].transpose(), extent = extentX, aspect = 'auto', cmap = cmap)
        _ = ax[2].imshow(XE['int'].transpose(), extent = extentY, aspect = 'auto', cmap = cmap)
        for a in ax: a.invert_yaxis()
        #
        ax[0].axvline(x = X - DX/2, linewidth = linewidth, color = 'red')
        ax[0].axvline(x = X + DX/2, linewidth = linewidth, color = 'red')
        ax[0].axhline(y = Y - DY/2, linewidth = linewidth, color = 'red')
        ax[0].axhline(y = Y + DY/2, linewidth = linewidth, color = 'red')
        for i in [1,2]: 
            ax[i].axhline(y = E - DE/2, linewidth = linewidth, color = 'red')
            ax[i].axhline(y = E + DE/2, linewidth = linewidth, color = 'red')
        #
        ax[0].set_xlabel('X (°)')
        ax[0].set_ylabel('Y (°)')
        ax[1].set_xlabel('Y (°)')
        ax[2].set_xlabel('X (°)')
        #
        fig.suptitle("ID {0}".format(D.get('Experiment', {}).get('Spectrum_ID', '')))

    Interact = ipw.interactive_output(plot, {'X': SliderX, 
                                             'Y': SliderY, 
                                             'E': SliderE, 
                                             'DX': SliderDX,
                                             'DY': SliderDY, 
                                             'DE': SliderDE})

    box_out = ipw.VBox([box_sliders, Interact])
    box_out.layout = ipw.Layout(border="solid 1px gray", margin="5px", padding="2")
    display(box_out)
    









# ==============================================================================================
# ==============================================================================================
# ==============================================================================================
# ==============================================================================================

def Explore(D = {}, **kwargs):
    """"""
    if not type(D) is dict: print(Fore.RED + 'Argument D must be a grumpy dict.' + Fore.RESET); return
    #
    recognized_kwargs = ['']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    dim = len(np.shape(D.get('int', np.array([]))))
    if dim < 2:
        print(Fore.RED + 'The data seem to be have too few dimensions. This method requires at least two axes.' + Fore.RESET); return
    elif dim == 2: 
        _Explore2d(D = D, **kwargs)
    elif dim == 3:
        _Explore3d(D = D, **kwargs)
    else:
        print(Fore.RED + 'The data seem to be have too many dimensions. This method handles maximum three axes.' + Fore.RESET); return


def _Explore2d(D = {}, **kwargs):
    """"""
    recognized_kwargs = ['figsize', 'cmap', 'rotate']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    figsize = kwargs.get('figsize', (5,3.5))
    cmap = kwargs.get('cmap', 'bone_r')
    rotate = kwargs.get('rotate', None)
    #
    # ------------------ get axes, labels, and types 
    axis1, axis1label, axis1type = D.get('x', np.array([])), D.get('Meta',{}).get('x_label', ''), D.get('Meta',{}).get('x_type', '')
    axis2, axis2label, axis2type = D.get('y', np.array([])), D.get('Meta',{}).get('y_label', ''), D.get('Meta',{}).get('y_type', '')
    #
    # ------------------- intensity
    intensity = D.get('int', np.array([]))
    # ------------------- rotate if angles are involved
    do_rotation = False
    if type(rotate) is type(None):
        if axis1type in [defax_angle, defax_xdefl, defax_ydefl] or axis2type in [defax_angle, defax_xdefl, defax_ydefl]:
            do_rotation = True
    else:
        do_rotation = rotate
    if do_rotation:
        axis1, axis2 = axis2, axis1
        axis1label, axis2label = axis2label, axis1label
        axis1type, axis2type = axis2type, axis1type
        intensity = intensity.transpose()
    # ------------------- intensity sliders
    c_min = intensity.min()
    c_max = intensity.max()
    c_step = (c_max-c_min)/100
    SliderMinInt = ipw.FloatSlider(min=c_min, max=c_max-c_step, step=c_step, description = 'Min.', value=c_min, readout_format = ".3f", continuous_update = True)
    SliderMaxInt = ipw.FloatSlider(min=c_min+c_step, max=c_max, step=c_step, description = 'Max.', value=c_max, readout_format = ".3f", continuous_update = True)
    def SliderMinInt_change(change):
        if SliderMinInt.value >= SliderMaxInt.value : SliderMaxInt.value = SliderMinInt.value + c_step
    SliderMinInt.observe(SliderMinInt_change, names='value')
    def SliderMaxInt_change(change):
        if SliderMaxInt.value <= SliderMinInt.value : SliderMinInt.value = SliderMaxInt.value - c_step
    SliderMaxInt.observe(SliderMaxInt_change, names='value')
    #
    # -------------------- axvline and axhline sliders
    x_step = (axis1.max()-axis1.min())/(len(axis1)-1)
    y_step = (axis2.max()-axis2.min())/(len(axis2)-1)
    SliderAxV = ipw.FloatSlider(min=axis1.min(), max=axis1.max(), step=x_step, description = 'X', value=axis1.min(), readout_format = ".3f", continuous_update = True)
    SliderAxH = ipw.FloatSlider(min=axis2.min(), max=axis2.max(), step=y_step, description = 'Y', value=axis2.min(), readout_format = ".3f", continuous_update = True)
    def SliderAxVH_change(change):
        pass
    SliderAxV.observe(SliderAxVH_change, names='value')
    SliderAxH.observe(SliderAxVH_change, names='value')
    #
    # ---------------------- the plot function
    def plot(vmin, vmax, axv, axh):
        fig, ax = plt.subplots(figsize = figsize)
        #
        extent = [axis1[0], axis1[-1], axis2[-1], axis2[0]]
        im = ax.imshow(intensity, extent = extent, vmin = vmin, vmax = vmax, cmap = cmap, aspect = 'auto')#aspect)
        if not axv == axis1.min():
            ax.axvline(x = axv)
        if not axh == axis2.min():
            ax.axhline(y = axh)
        ax.invert_yaxis()
        ax.set_xlabel(axis1label)
        ax.set_ylabel(axis2label)
        ax.set_title("ID {0}".format(D.get('Experiment', {}).get('Spectrum_ID', '')))
        plt.tight_layout()
    #
    # --------------------- interctive stuff
    Interact = ipw.interactive_output(plot, {'vmin': SliderMinInt, 
                                             'vmax': SliderMaxInt,
                                             'axv': SliderAxV,
                                             'axh': SliderAxH})
    box_1 = ipw.HBox([ipw.VBox([SliderMinInt, SliderMaxInt]), ipw.VBox([SliderAxV, SliderAxH])])
    box = ipw.VBox([Interact, box_1])
    display(box)


def _Explore3d(D = {}, **kwargs):
    """"""
    recognized_kwargs = ['figsize', 'cmap']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    figsize = kwargs.get('figsize', (4,3))
    cmap = kwargs.get('cmap', 'bone_r')
    # ----------------- get axes and stuff
    xaxis = D.get('x')
    xlabel = D.get('Meta',{}).get('x_label')
    xtype = D.get('Meta',{}).get('x_type')
    yaxis = D.get('y')
    ylabel = D.get('Meta',{}).get('y_label')
    ytype = D.get('Meta',{}).get('y_type')
    zaxis = D.get('z')
    zlabel = D.get('Meta',{}).get('z_label')
    ztype = D.get('Meta',{}).get('z_type')
    intensity = D.get('int')
    #
    # ----------------- start conditions
    sum_axis = xaxis; sum_axis_label = xlabel; sum_axis_type = xtype
    axis1 = yaxis;    axis1_label =    ylabel; axis1_type =    ytype
    axis2 = zaxis;    axis2_label =    zlabel; axis2_type =    ztype
    #
    # ----------------- method to get the non-sum_axis axes
    def getaxs(sumaxislbl):
        if sumaxislbl == xlabel:
            return D.get('y'), D.get('Meta',{}).get('y_label'), D.get('z'), D.get('Meta',{}).get('z_label')
        elif sumaxislbl == ylabel:
            return D.get('x'), D.get('Meta',{}).get('x_label'), D.get('z'), D.get('Meta',{}).get('z_label')
        elif sumaxislbl == zlabel:
            return D.get('x'), D.get('Meta',{}).get('x_label'), D.get('y'), D.get('Meta',{}).get('y_label')
    #
    # ----------------- method to make a 2d slice
    def makeslice(axis, v, dv):
        v1 = v - dv/2
        v2 = v + dv/2
        if axis == xlabel: 
            v1, v2, iv1, iv2 = _fixIndices(xaxis, v1, v2, None, None)
            return np.array(intensity[:,:,iv1:iv2+1].sum(axis = 2))
        elif axis == ylabel:
            v1, v2, iv1, iv2 = _fixIndices(yaxis, v1, v2, None, None)
            return np.array(intensity[:,iv1:iv2+1,:].sum(axis = 1))
        elif axis == zlabel:
            v1, v2, iv1, iv2 = _fixIndices(zaxis, v1, v2, None, None)
            return np.array(intensity[iv1:iv2+1,:,:].sum(axis = 0))
    #
    # ----------------- intensity slide bars
    c_min = intensity.min()
    c_max = intensity.max()
    c_step = (c_max-c_min)/100
    SliderMinInt = ipw.FloatSlider(min=c_min, max=c_max-c_step, step=c_step, description = 'Min.', value=c_min, readout_format = ".3f", continuous_update = True)
    SliderMaxInt = ipw.FloatSlider(min=c_min+c_step, max=c_max, step=c_step, description = 'Max.', value=c_max, readout_format = ".3f", continuous_update = True)
    def SliderMinInt_change(change):
        if SliderMinInt.value >= SliderMaxInt.value : SliderMaxInt.value = SliderMinInt.value + c_step
    SliderMinInt.observe(SliderMinInt_change, names='value')
    def SliderMaxInt_change(change):
        if SliderMaxInt.value <= SliderMinInt.value : SliderMinInt.value = SliderMaxInt.value - c_step
    SliderMaxInt.observe(SliderMaxInt_change, names='value')
    x_step = (axis1.max()-axis1.min())/(len(axis1)-1)
    y_step = (axis2.max()-axis2.min())/(len(axis2)-1)
    SliderAxV = ipw.FloatSlider(min=axis1.min(), max=axis1.max(), step=x_step, description = 'X', value=axis1.min(), readout_format = ".3f", continuous_update = True)
    SliderAxH = ipw.FloatSlider(min=axis2.min(), max=axis2.max(), step=y_step, description = 'Y', value=axis2.min(), readout_format = ".3f", continuous_update = True)
    #
    # ----------------- select sum_axis
    RadioButtons_Axis = ipw.RadioButtons(
                                options=[xlabel, ylabel, zlabel],
                                value=xlabel,
                                description='Axes',
                                disabled=False
                                )
    def RadioButtons_Axis_change(change):
        if RadioButtons_Axis.value == xlabel:
            sum_axis = xaxis
        elif RadioButtons_Axis.value == ylabel:
            sum_axis = yaxis
        else: 
            sum_axis = zaxis
        axis1, axis1_label, axis2, axis2_label = getaxs(RadioButtons_Axis.value)
        if sum_axis.max() >=  SliderAxisValue.min:
            SliderAxisValue.max = sum_axis.max()
            SliderAxisValue.min = sum_axis.min()
        else:
            SliderAxisValue.min = sum_axis.min()
            SliderAxisValue.max = sum_axis.max()
        SliderAxisValue.value = sum_axis.mean()
        intensity = makeslice(RadioButtons_Axis.value, sum_axis.mean(), (sum_axis.max()-sum_axis.min())/(len(sum_axis)-1))
        c_min = intensity.min()
        c_max = intensity.max()
        c_step = (c_max-c_min)/100.
        SliderMinInt.min, SliderMinInt.max, SliderMinInt.step, SliderMinInt.value = c_min, c_max - c_step, c_step, c_min
        SliderMaxInt.min, SliderMaxInt.max, SliderMaxInt.step, SliderMaxInt.value = c_min + c_step, c_max, c_step, c_max
        if axis1.max() > SliderAxV.min:
            SliderAxV.max = axis1.max()
            SliderAxV.min = axis1.min()
        else:
            SliderAxV.min = axis1.min()
            SliderAxV.max = axis1.max()
        if axis2.max() > SliderAxH.min:
            SliderAxH.max = axis2.max()
            SliderAxH.min = axis2.min()
        else:
            SliderAxH.min = axis2.min()
            SliderAxH.max = axis2.max()
        SliderAxV.step, SliderAxV.value = (axis1.max()-axis1.min())/(len(axis1)-1), axis1.min()
        SliderAxH.step, SliderAxH.value = (axis2.max()-axis2.min())/(len(axis2)-1), axis2.min()
    RadioButtons_Axis.observe(RadioButtons_Axis_change, names='value')
    #
    # ----------------- sum_axis value sliders
    SliderAxisValue = ipw.FloatSlider(min=sum_axis.min(), 
                                      max=sum_axis.max(), 
                                      step=(sum_axis.max()-sum_axis.min())/(len(sum_axis)-1), 
                                      description = 'Value', 
                                      value=sum_axis.mean(), 
                                      readout_format = ".3f", 
                                      continuous_update = True)
    SliderAxisDelta = ipw.FloatSlider(min=(sum_axis.max()-sum_axis.min())/(len(sum_axis)-1), 
                                      max=sum_axis.max()-sum_axis.min(), 
                                      step=(sum_axis.max()-sum_axis.min())/(len(sum_axis)-1), 
                                      description = 'Delta', 
                                      value=(sum_axis.max()-sum_axis.min())/(len(sum_axis)-1), 
                                      readout_format = ".3f", 
                                      continuous_update = True)
    def SliderAxis_change(change):
        global intensity
        #intensity = slice(RadioButtons_Axis.value, SliderAxisValue.value, SliderAxisDelta.value)
    SliderAxisValue.observe(SliderAxis_change, names='value')
    SliderAxisDelta.observe(SliderAxis_change, names='value')
    #
    # ----------------- plot
    def plot(vmin, vmax, axis_value, axis_delta, sum_axis, axv, axh):
        fig, ax = plt.subplots(figsize = figsize)

        axis1, axis1label, axis2, axis2label = getaxs(sum_axis)
        slice = makeslice(sum_axis, axis_value, axis_delta)
        
        extent = [axis1[0], axis1[-1], axis2[-1], axis2[0]]
        im = ax.imshow(slice, extent = extent, vmin = vmin, vmax = vmax, cmap = cmap, aspect = 'auto')
        ax.invert_yaxis()
        if not axv == axis1.min():
            ax.axvline(x = axv)
        if not axh == axis2.min():
            ax.axhline(y = axh)
        ax.set_xlabel(axis1label)
        ax.set_ylabel(axis2label)
        ax.set_title("ID {0}".format(D.get('Experiment', {}).get('Spectrum_ID', '')))
        plt.tight_layout()
    #
    # ----------------- interact stuff
    Interact = ipw.interactive_output(plot, {'vmin': SliderMinInt, 
                                             'vmax': SliderMaxInt,
                                             'axis_value': SliderAxisValue,
                                             'axis_delta': SliderAxisDelta,
                                             'sum_axis': RadioButtons_Axis,
                                             'axv': SliderAxV,
                                             'axh': SliderAxH})
    #
    box_sliders_int = ipw.VBox([SliderMinInt, SliderMaxInt])
    box_sliders_axvh = ipw.VBox([SliderAxV, SliderAxH])
    box_sliders1 = ipw.HBox([box_sliders_int, box_sliders_axvh])
    #
    box_1 = ipw.VBox([Interact, box_sliders1])
    box_2 = ipw.VBox([RadioButtons_Axis, SliderAxisValue, SliderAxisDelta])
    box_out = ipw.HBox([box_1, box_2])
    display(box_out)



        




            


# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================

class SpinEDC():
    """
    Interactive manipulation of spin edc data.
    Use:
        interactive = SpinEDC(loaded_dict)
        new_dict = interactive.result()
    """
    def __init__(self, data = {}, s = None):
        if not type(data) is dict:
            print(Fore.RED + 'Argument data must a grumpy dict.' + Fore.RESET); return
        if not data.get('Type', '') == 'spin_edc':
            print(Fore.RED + 'The dicts does not contain spin edc data.' + Fore.RESET); return
        #
        print(Fore.BLUE + "Deselect edcs you don't want to use. Normalize the edcs at a certain energy.")
        print("Export the result (to a new dict) using the SpinEDC.result() method.\n" + Fore.RESET)
        #
        if type(s) is type(None):
            global SHERMAN
            self.sherman = SHERMAN
        else:
            self.sherman = s
        self.data = data
        self.Data = data.get("Data")
        self.Parameter_values = data.get("Parameter_values")[0]
        self.energy = self.data.get("x", np.array([]))

        self.all_curves = []
        self.avg_on  = []
        self.avg_off = []
        self.asym = []
        self.included = []

        self.include_checkboxes = []
        self.interact_args = {}
        for i, par in enumerate(self.Parameter_values):
            if "OFF" in par: mpol = f"off {i+1}"
            else: mpol = f"on {i+1}"
            self.include_checkboxes.append(ipw.Checkbox(value = True, description = mpol, width = 10))
            self.interact_args.update({mpol: self.include_checkboxes[-1]})
        self._box1 = ipw.HBox(layout = ipw.Layout(width="100%", display="inline-flex", flex_flow="row_wrap"))
        self._box1.children = [chbx for chbx in self.include_checkboxes]
        #self._box1 = ipw.HBox(self.include_checkboxes)
        
        self.legend1_checkbox = ipw.Checkbox(value = True, description = "Show legend")
        self.interact_args.update({"legend1": self.legend1_checkbox})
        self._box2 = ipw.HBox([self.legend1_checkbox])

        self.norm_checkbox = ipw.Checkbox(value = False, description = "Normalize")
        self.norm_energy_slider = ipw.FloatSlider(min = self.Data[0][0][0], max = self.Data[0][0][-1], 
            step = self.Data[0][0][1] - self.Data[0][0][0], description = 'Energy', readout_format = ".000f")
        self.norm_energy_width_slider = ipw.FloatSlider(min = self.Data[0][0][1] - self.Data[0][0][0], 
            max = 10 * (self.Data[0][0][1] - self.Data[0][0][0]), step = self.Data[0][0][1] - self.Data[0][0][0], description = 'Width')
        self.interact_args.update({"normalize": self.norm_checkbox, "norm_energy": self.norm_energy_slider, "norm_energy_width": self.norm_energy_width_slider})
        self._box3 = ipw.VBox([self.norm_checkbox, self.norm_energy_slider, self.norm_energy_width_slider])

        self.sherman_slider = ipw.FloatSlider(min = 0.05, max = 0.65, step = 0.01, description = 'Sherman', value = self.sherman)
        self.interact_args.update({"sherman": self.sherman_slider})
        self._box4 = ipw.HBox([self._box3, self.sherman_slider])

        def Plot(**kwargs):
            fig, ax = plt.subplots(ncols = 3, figsize = (12,4))
            #
            self.Parameter_values_used = []
            self.all_curves = []
            for d in self.Data: self.all_curves.append(d[1])
            #
            if self.norm_checkbox.value:
                ax[0].axvline(x = self.norm_energy_slider.value, color = "k", linestyle = "--", linewidth = 0.65)
                ax[0].axvline(x = self.norm_energy_slider.value + self.norm_energy_width_slider.value, color = "k", linestyle = "--", linewidth = 0.65)
                i1 = abs(self.norm_energy_slider.value - self.Data[0][0]).argmin()
                i2 = abs(self.norm_energy_slider.value + self.norm_energy_width_slider.value - self.Data[0][0]).argmin()+1
                for i, d in enumerate(self.all_curves):
                    norm = d[i1:i2].sum()/(i2-i1)
                    self.all_curves[i] = d/norm
            #
            self.avg_on, self.avg_off =  np.zeros(len(self.Data[0][0])), np.zeros(len(self.Data[0][0]))
            n_on, n_off = 0, 0
            self.included = []
            for i, d in enumerate(self.all_curves):
                if "OFF" in self.Parameter_values[i]:
                    mpol = f"off {i+1}"
                    if self.include_checkboxes[i].value:
                        self.avg_off += d
                        self.included.append(i)
                        n_off += 1
                        liwi = 1; listy = "-"
                    else:
                        liwi = 0.75; listy = ":"
                else:
                    mpol = f"on {i+1}"
                    if self.include_checkboxes[i].value:
                        self.avg_on += d
                        self.included.append(i)
                        n_on += 1
                        liwi = 1; listy = "-"
                    else:
                        liwi = 0.75; listy = ":"
                ax[0].plot(self.energy, d, linestyle = listy, linewidth = liwi, label = mpol)
            #
            if n_off > 0: self.avg_off /= n_off
            if n_on > 0: self.avg_on /= n_on
            ax[1].plot(self.energy, self.avg_off, color = "blue", label = "off")
            ax[1].plot(self.energy, self.avg_on,  color = "red",  label = "on")
            #
            self.asym = np.zeros(len(self.energy))
            for i in range(len(self.asym)):
                denom = self.avg_off[i] + self.avg_on[i]
                if not denom < 1e-4:
                    self.asym[i] = (self.avg_off[i] - self.avg_on[i])/denom
                else:
                    self.asym[i] = 0
            self.sherman = self.sherman_slider.value
            self.asym /= self.sherman
            ax[2].plot(self.energy, self.asym, color = 'k')
            for i, txt in enumerate(["EDCs", "Average", "Asymmetry"]):
                ax[i].set_title(txt)
                ax[i].set_xlabel("Energy, eV")
            if self.legend1_checkbox.value:
                ax[0].legend()
                ax[1].legend()
            fig.tight_layout()

        self._Interact = ipw.interactive_output(Plot, self.interact_args)
        self._box_out = ipw.VBox([self._box1, self._box2, self._Interact, self._box4])
        display(self._box_out)

    def result(self):
        result = deepcopy(self.data)
        #
        pvs = []
        for i, pv in enumerate(self.Parameter_values):
            if i in self.included: pvs.append(pv)
        result.update({"Parameter_values": pvs})
        #
        Data = []
        for i, edc in enumerate(self.all_curves):
            if i in self.included:
                Data.append(np.array([self.energy, edc]))
        result.update({"Data": Data})
        #
        result.update({"int": self.avg_off, "int_on": self.avg_on, "asymmetry": self.asym})
        #
        all_off, all_on = [], []
        for i, edc in enumerate(self.all_curves):
            if i in self.included:
                if "OFF" in self.Parameter_values[i]:
                    all_off.append(edc)
                else:
                    all_on.append(edc)
        result.update({"int_all": all_off, "int_all_on": all_on})
        #
        return result





 # ============================================================================================================   


def QuickSpin(D = {}, **kwargs):

    """
    For manipulation of simple spin edcs.

    Args:
        sherman (float)
        exclude (list, indices of scans to exclude)
        median_filter (bool), size (int), mode (str)
        remove_spikes (bool), threshold (float)
            filter_outliers (LEGACY. Same as remove_spikes. Will be removed)
        plot (bool), figsize (tupple), linewidth (float), linestyle (str)
        shup (bool)
        norm (bool, together with norm_pnts)

    """
    #
    recognized_kwargs = ['plot', 'shup', 'figsize', 'linewidth', 'linestyle', 'exclude', 'median_filter', 
                         'size', 'mode', 'remove_spikes', 'filter_outliers', 'threshold', 'sherman',
                         'norm', 'norm_start_point', 'norm_points']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    plot = kwargs.get('plot', True)
    figsize = kwargs.get('figsize', (10,2.5))
    linewidth = kwargs.get('linewidth', 0.75)
    linestyle = kwargs.get('linestyle', '-')
                           
    exclude = kwargs.get('exclude', [])
    medianfilter = kwargs.get('median_filter', False)
    size = kwargs.get('size', 3)
    mode = kwargs.get('mode', 'reflect')
    remove_spikes = kwargs.get('remove_spikes', None)
    if type(remove_spikes) is type(None):
        remove_spikes = kwargs.get('filter_outliers', False)        
    threshold = kwargs.get('threshold', 2)

    norm = kwargs.get("norm", False)
    norm_start_point = int(kwargs.get("norm_start_point", 0))
    norm_points = (kwargs.get("norm_points", 5))
    if norm:
        if norm_start_point < 0: norm_start_point = 0
        if norm_points <1: norm_points = 1  

    global SHERMAN
    sherman = kwargs.get('sherman', SHERMAN)

    shup = kwargs.get('shup', False)
    #
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must a grumpy dict.' + Fore.RESET); return {}
    if not D.get('Type', '') == 'spin_edc':
        print(Fore.RED + 'The dicts does not contain spin edc data.' + Fore.RESET); return D

    EDC_OFF = D.get('int_all')
    EDC_ON = D.get('int_all_on')

    if plot:
        fig, ax = plt.subplots(figsize = figsize, ncols = 4)
        plt.tight_layout()
        for a in ax: a.set_xlabel('Energy (ev)')

    energy = D.get('x')
    edc_off = np.copy(energy) * 0
    edc_on = np.copy(edc_off)

    i, ion, ioff = 0, 0, 0
    for edc in EDC_OFF:
        if not i in exclude:
            if plot: ax[0].plot(energy, edc, linewidth = linewidth, color = 'tab:blue', linestyle = linestyle)
            edc_off += edc
            ioff += 1
        i += 1
    for edc in EDC_ON:
        if not i in exclude:
            if plot: ax[0].plot(energy, edc, linewidth = linewidth, color = 'tab:red', linestyle = linestyle)
            edc_on += edc
            ion += 1
        i += 1
    edc_off /= ioff
    edc_on /= ion

    if remove_spikes:
        edc_off = removeSpikes(edc_off, threshold = threshold)
        edc_on =  removeSpikes(edc_on, threshold = threshold)
        if not shup: print('Applied removeSpikes()')

    if medianfilter:
        edc_off = medianFilter(edc_off, size = size, mode = mode)
        edc_on =  medianFilter(edc_on, size = size, mode = mode)
        if not shup: print('Applied medianFilter()')
    
    if norm:
        edc_on = edc_on / edc_on[norm_start_point:norm_start_point + norm_points].sum() * edc_off[norm_start_point:norm_start_point + norm_points].sum()
    
    if plot:
        ax[1].plot(energy, edc_off, linewidth = linewidth, color = 'tab:blue', label = 'pol. +', linestyle = linestyle)
        ax[1].plot(energy, edc_on, linewidth = linewidth, color = 'tab:red', label = 'pol. -', linestyle = linestyle)
        ax[1].legend()
        ax[0].set_title('edc')
        ax[1].set_title('edc avg.')
    
    asymmetry = np.zeros([len(energy)])
    for i in range(len(energy)):
        denom = (edc_off[i] + edc_on[i])
        if not denom == 0:
            nom = (edc_off[i] - edc_on[i])
            asymmetry[i] = (nom/denom)
        else:
            asymmetry[i] = 0
    if plot:
        ax[2].plot(energy, asymmetry / sherman, linewidth = linewidth, color = 'k')
        ax[2].set_title(f'asymmetry (S = {sherman})')
    
    P1 = (edc_off + edc_on) * (1 + asymmetry / sherman) * 0.5 # there must be a 0.5 here, right?
    P2 = (edc_off + edc_on) * (1 - asymmetry / sherman) * 0.5
    if plot:
        ax[3].plot(energy, P1, color = 'tab:blue', label = 'P1', linewidth = linewidth, linestyle = linestyle)
        ax[3].plot(energy, P2, color = 'tab:red', label = 'P2', linewidth = linewidth, linestyle = linestyle)
        ax[3].legend()
        ax[3].set_title('Components, S = {0}'.format(sherman))
    
    DD = copy.deepcopy(D)
    DD.update({'int': edc_off})
    DD.update({'int_on': edc_on})
    DD.update({'asymmetry': asymmetry/sherman})
    return DD



# ==============================================================================================


def QuickSpinMDC(D = {}, **kwargs):
    """
    For manipulation of deflector spin data.

    kwargs:
        sherman (float)
        median_filter (bool), size (int), mode (str)
        remove_spikes (bool), threshold (float)
            filter_outliers (LEGACY. Same as remove_spikes. Will be removed)
        plot (bool), figsize (tupple), linewidth (float), linestyle (str)
            vmin (float), vmax (float), vmina (float), vmaxa (float), cmap (str), colorbar (bool)
        shup (bool)
        
    """
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must a grumpy dict.' + Fore.RESET); return D
    if not D.get('Type', '').startswith('spin_deflector'):
        print(Fore.RED + 'The dicts does not contain deflector spin data.' + Fore.RESET); return D
    #
    recognized_kwargs = ['plot', 'shup', 'median_filter_edc', 'size', 'mode', 'remove_spikes', 'filter_outliers',
                        'threshold', 'sherman', 'figsize', 'linewidth', 'linestyle', 'vmin', 'vmax', 'vmina',
                        'vmaxa', 'cmap', 'colorbar', 'median_filter']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    plot = kwargs.get('plot', True)
    shup = kwargs.get('shup', False)
    medianfilter = kwargs.get('median_filter_edc', False)
    size = kwargs.get('size', 3)
    mode = kwargs.get('mode', 'reflect')
    remove_spikes = kwargs.get('remove_spikes', None)
    if type(remove_spikes) is type(None):
        remove_spikes = kwargs.get('filter_outliers', False)
    threshold = kwargs.get('threshold', 2)
    
    global SHERMAN
    sherman = kwargs.get('sherman', SHERMAN)
    figsize = kwargs.get('figsize', (0,0))
    linewidth = kwargs.get('linewidth', 0.75)
    linestyle = kwargs.get('linestyle', '-')
    vmin = kwargs.get('vmin', None)
    vmax = kwargs.get('vmax', None)
    vmina = kwargs.get('vmina', None)
    vmaxa = kwargs.get('vmaxa', None)
    cmap = kwargs.get('cmap', 'bone_r')
    colorbar = kwargs.get('colorbar', False)


    # Note: change the conditions below when incorporating 2d deflection scans...
    energy  = D.get('x')
    defl_x  = D.get('y', np.array([]))
    defl_y  = D.get('z', np.array([]))
    edc_off = D.get('int', np.array([]))
    edc_on  = D.get('int_on', np.array([]))
    asym    = D.get('asymmetry', np.array([]))
    if len(defl_x) > 0 and len(defl_y) == 0:
        direction = 'x'
        defl_axis = np.copy(defl_x)
    elif len(defl_x) == 0 and len(defl_y) > 0: 
        direction = 'y'
        defl_axis = np.copy(defl_y)
    elif len(defl_x) > 0 and len(defl_y) > 0:
        direction = 'xy'
        defl_axis = np.array(range(len(edc_off)))
    else:
        print(Fore.RED + 'Probably most likely there is potentially some sloppy coding involved here. Sorry. Aborting.' + Fore.RESET); return {}        

    # energy is an array with one vaule if fixed energy, otherwise several values
    if len(energy) == 1: EDC = False
    else: EDC = True

    #
    if not EDC:
        if removeSpikes:
            edc_off = removeSpikes(edc_off, threshold = threshold)
            edc_on =  removeSpikes(edc_on, threshold = threshold)
            if not shup: print('Applied removeSpikes()')
        if medianfilter:
            edc_off = medianFilter(edc_off, size = size, mode = mode)
            edc_on =  medianFilter(edc_on, size = size, mode = mode)
            if not shup: print('Applied medianFilter()')
        #
        asymmetry = np.zeros([len(edc_off)])
        for i in range(len(edc_off)):
            denom = (edc_off[i] + edc_on[i])
            if not denom == 0:
                nom = (edc_off[i] - edc_on[i])
                asymmetry[i] = (nom/denom)
            else:
                asymmetry[i] = 0
        #
        P1 = (edc_off + edc_on) * (1 + asymmetry / sherman) * 0.5 # there must be a 0.5 here, right?
        P2 = (edc_off + edc_on) * (1 - asymmetry / sherman) * 0.5
        #
        if plot:
            if figsize == (0,0): figsize = (10,3)
            fig, ax = plt.subplots(figsize = figsize, ncols = 3)
            ax[0].plot(defl_axis, edc_off, color = 'tab:blue', linewidth = linewidth, label = 'Off', linestyle = linestyle)
            ax[0].plot(defl_axis, edc_on, color = 'tab:red', linewidth = linewidth, label = 'On', linestyle = linestyle)
            ax[1].plot(defl_axis, asymmetry, color = 'k', linewidth = linewidth, label = 'Asymmetry', linestyle = linestyle)
            ax[2].plot(defl_axis, P1, color = 'tab:blue', linewidth = linewidth, label = 'P1', linestyle = linestyle)
            ax[2].plot(defl_axis, P2, color = 'tab:red', linewidth = linewidth, label = 'P2', linestyle = linestyle)
            #
            if direction in ['x','y']:
                xlabel = 'Deflection {0} (deg.)'.format(direction)
            else:
                xlabel = 'Defl. ({0},{1}) to ({2},{3}) (deg.)'.format(defl_x[0], defl_y[0], defl_x[-1], defl_y[-1])
            #
            for i, ttl in enumerate(['Intensity', 'Asymmetry', 'Intensity (S={0})'.format(sherman)]):
                ax[i].set_title(ttl)
                ax[i].set_xlabel(xlabel)
                if i in [0,2]:
                    ax[i].legend()
                    ax[i].set_ylabel('Intensity')
                if i == 1: ax[i].set_ylabel('Asymmetry')
            plt.tight_layout()
        DD = copy.deepcopy(D)
        DD.update({'int': edc_off})
        DD.update({'int_on': edc_on})
        DD.update({'asymmetry': asymmetry})
        
    else:
        off_map, on_map, asym_map, p1_map, p2_map = [], [], [], [], []
        for i in range(len(asym)):
            edc_off_f = edc_off[i]
            edc_on_f  = edc_on[i]
            if removeSpikes:
                edc_off_f = removeSpikes(edc_off_f, threshold = threshold)
                edc_on_f  = removeSpikes(edc_on_f, threshold = threshold)
            if medianfilter:
                edc_off_f = medianFilter(edc_off_f, size = size, mode = mode)
                edc_on_f =  medianFilter(edc_on_f, size = size, mode = mode)
            asymmetry = (edc_off_f - edc_on_f) / (edc_off_f + edc_on_f)
            p1 = (edc_off_f + edc_on_f) * (1 + asymmetry / sherman) * 0.5 # there must be a 0.5 here, right?
            p2 = (edc_off_f + edc_on_f) * (1 - asymmetry / sherman) * 0.5
            off_map.append(edc_off_f)
            on_map.append(edc_on_f)
            asym_map.append(asymmetry)
            p1_map.append(p1)
            p2_map.append(p2)
        off_map, on_map, asym_map, p1_map, p2_map = np.array(off_map), np.array(on_map), np.array(asym_map), np.array(p1_map), np.array(p2_map)
        #
        if plot:
            if figsize == (0,0): figsize = (11,7)
            fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize = figsize)
            ax = ax.flatten()
            ax[5].remove()
            extent = [energy[0], energy[-1], defl_axis[-1], defl_axis[0]]
            IMS = []
            IMS.append( ax[0].imshow(off_map, extent = extent, aspect = 'auto', vmin = vmin, vmax = vmax, cmap = cmap) )
            IMS.append( ax[1].imshow(on_map, extent = extent, aspect = 'auto', vmin = vmin, vmax = vmax, cmap = cmap) )
            IMS.append( ax[2].imshow(asym_map, extent = extent, aspect = 'auto', vmin = vmina, vmax = vmaxa, cmap = cmap) )
            IMS.append( ax[3].imshow(p1_map, extent = extent, aspect = 'auto', vmin = vmin, vmax = vmax, cmap = cmap) )
            IMS.append( ax[4].imshow(p2_map, extent = extent, aspect = 'auto', vmin = vmin, vmax = vmax, cmap = cmap) )
            ttls = ['Intensity Off', 'Intensity On', 'Asymmetry', 'Intensity S={0}'.format(sherman), 'Intensity S={0}'.format(sherman)]
            for i, ttl in enumerate(ttls):
                ax[i].set_title(ttl) 
                if colorbar: fig.colorbar(IMS[i], ax = ax[i])
                if i in [2,3,4]: ax[i].set_xlabel('Energy (eV)')
                if i in [0,3]: ax[i].set_ylabel('Deflection')
                ax[i].invert_yaxis()
        DD = copy.deepcopy(D)
        DD.update({'int': off_map})
        DD.update({'int_on': on_map})
        DD.update({'asymmetry': asym_map})
    
    return DD

        







# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================

# not tested after re-write

def Polarization(D1 = {'Measurement_type': 'none'}, D2 = {'Measurement_type': 'none'}, D3 = {'Measurement_type': 'none'}, **kwargs):
    """
    Calculates Px, Py, and Pz from three measurements, Px and Py from two measurements, or Pz from one,
    given that the measurements are from the correct coil and rotor combinations.

    Arguments:
        D1      Coil 2, Rotor +1
        D2      Coil 2, Rotor -1
        D3      Coil 1, Rotor +/-1
    """
    #print('This method is written but not tested for all variations of data. Report errors. Thanks.\n')
    #
    recognized_kwargs = ['SF', 'plot', 'figsize', 'cmap', 'vmin', 'vmax']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)
    #
    global SHERMAN
    SF = kwargs.get('sherman', SHERMAN)
    plot = kwargs.get('plot', True)
    figsize = kwargs.get('figsize', (4,3))
    cmap = kwargs.get('cmap', 'rainbow')
    vmin = kwargs.get('vmin', None)
    vmax = kwargs.get('vmax', None)
    #
    allowed_types = ['spin_edc', 'spin_deflector_1a', 'spin_deflector_1b']
    #
    D = [D1, D2, D3]
    for i,d in enumerate(D):
        if not type(d) is dict: D[i] = {}
    #
    types = []
    for d in D: types.append(d.get('Type', 'none'))
    #
    if all(x=='none' for x in types):
        print(Fore.GREEN + 'Please use help(grumpy.Polarization).' + Fore.RESET)
        return {}
    for i,t in enumerate(types):
        if not t == 'none' and not t in allowed_types: types[i] = 'wrong'
    if any(x=='wrong' for x in types):
        print(Fore.RED + 'One or more of the passed data dicts are of the wrong measurement types for this method.' + Fore.RESET)
        return {}
    #
    Case = ''
    if all(types[0] == t for t in types) and not types[0] == 'none':    # Px, Py, Pz
        Case = 'Pxyz'
    elif types[0] == types[1] and not types[0] == 'none':               # Px, Py
        Case = 'Pxy'
    elif not types[2] == 'none':                                        # Pz
        Case = 'Pz'
    else:
        print(Fore.GREEN + 'For just one data set with coil 2 use QuickSpin() or QuickSpinMDC().' + Fore.RESET)
    #
    if Case in ['Pxyz', 'Pxy']: j = 0
    else: j = 2
    #
    
    #
    results = {}
    Px, Py, Pz = None, None, None
    # SPIN_EDC (w/o deflectors)
    if any(x=='spin_edc' for x in types):
        if plot: fig, ax = plt.subplots(figsize = figsize)
        if Case == 'Pxyz' or Case == 'Pxy':
            Px, Py = _polarization_xy(D[0]['asymmetry'], D[1]['asymmetry'], SF)
        elif Case == 'Pxyz' or Case == 'Pz':
            Pz = _polarization_z(D[2]['asymmetry'], SF)
        else:
            pass
        results.update({'energy': D[0].get('x')})
        if not type(Px) is type(None):
            results.update({'px': Px})
            if plot: ax.plot(results['energy'], Px, label = 'Px')
        if not type(Py) is type(None):
            results.update({'py': Py})
            if plot: ax.plot(results['energy'], Py, label = 'Py')
        if not type(Pz) is type(None):
            results.update({'pz': Pz})
            if plot: ax.plot(results['energy'], Pz, label = 'Pz')
        if plot:
            ax.set_xlabel(D[0].get('Meta',{}).get('x_label'))
            ax.set_ylabel('Polarization')
            ax.legend()
    # spin with deflectors
    else:
        # sort out which deflectors were used
        deflector1 = D.get('y')
        if D.get('Meta',{}).get('y_label','').lower().startswith('x'):
            dlbl1 = 'X'
            results.update({'x_deflection': deflector1})
        else: 
            dlbl1 = 'Y'
            results.update({'x_deflection': deflector1})
        xlabel = "{0}-deflection (deg.)".format(dlbl1)
        deflector = np.copy(deflector1)
        if 'y2' in D:
            deflector2 = D.get('y2')
            if D.get('Meta',{}).get('y2_label','').lower().startswith('x'):
                dlbl2 = 'X'
                results.update({'x_deflection': deflector2})
            else: 
                dlbl2 = 'Y'
                results.update({'y_deflection': deflector2})
            xlabel = "{0}{1}-deflection ({2},{3}) - ({4},{5}) (deg.)".format(dlbl1, dlbl2, deflector1[0], deflector2[0], deflector1[-1], deflector2[-1])
            deflector = np.array(range(len(deflector2))) # or len(deflector1), same length
        # Fixed energy
        if len(D[0]['x']) == 1:
            if plot: fig, ax = plt.subplots(figsize = figsize)
            if Case == 'Pxyz' or Case == 'Pxy':
                Px, Py = _polarization_xy(D[0]['asymmetry'], D[1]['asymmetry'], SF)
                results.update({'px': Px})
                results.update({'py': Py})
            elif Case == 'Pxyz' or Case == 'Pz':
                Pz = _polarization_z(D[2]['asymmetry'], SF)
                results.update({'pz': Pz})
            else:
                pass # should not happen...
            if not type(Px) is type(None):
                if plot: ax.plot(deflector, Px, label = 'Px')
            if not type(Py) is type(None):
                if plot: ax.plot(deflector, Py, label = 'Py')
            if not type(Pz) is type(None):
                if plot: ax.plot(deflector, Pz, label = 'Pz')
            if plot:
                ax.set_xlabel(xlabel)
                ax.set_ylabel('Polarization')
                ax.legend()
        # EDCs
        else:
            results.update({'energy': np.array(D[0]['x'])})
            pxmap, pymap, pzmap = [], [], []
            if Case == 'Pxyz' or Case == 'Pxy':
                for i in range(len(D[0]['asymmetry'])):
                    a0 = D[0]['asymmetry'][i]
                    a1 = D[1]['asymmetry'][i]
                    px, py = _polarization_xy(a0, a1, SF)
                    pxmap.append(px)
                    pymap.append(py)
                pxmap = np.array(pxmap)
                pymap = np.array(pymap)
                results.update({'px': pxmap})
                results.update({'py': pymap})
            elif Case == 'Pxyz' or Case == 'Pz':
                for i in range(len(D[2]['asymmetry'])):
                    a2 = D[2]['asymmetry'][i]
                    pz = _polarization_z(a2, SF)
                    pzmap.append(pz)
                pzmap = np.array(pzmap)
                results.update({'pz': pzmap})
            if plot:
                extent = [results['energy'][0], results['energy'][-1], deflector[-1], deflector[0]]
                if not type(pxmap) is list and not type(pymap) is list and not type(pzmap) is list:
                    fig, ax = plt.subplots(figsize = (figsize[0]*3, figsize[1]), ncols = 3)
                    im0 = ax[0].imshow(pxmap, extent = extent, aspect = 'auto', cmap = cmap, vmin = vmin, vmax = vmax)
                    im1 = ax[1].imshow(pymap, extent = extent, aspect = 'auto', cmap = cmap, vmin = vmin, vmax = vmax)
                    im2 = ax[2].imshow(pzmap, extent = extent, aspect = 'auto', cmap = cmap, vmin = vmin, vmax = vmax)
                    for i, t in enumerate(['Px', 'Py', 'Pz']):
                        ax[i].invert_yscale()
                        ax[i].set_title(t)
                        ax[i].set_xlabel(D[0].get('Meta',{}).get('x_label', ''))
                        ax[i].set_ylabel(xlabel)
                elif not type(pxmap) is list and not type(pymap) is list:
                    fig, ax = plt.subplots(figsize = (figsize[0]*2, figsize[1]), ncols = 2)
                    im0 = ax[0].imshow(pxmap, extent = extent, aspect = 'auto', cmap = cmap, vmin = vmin, vmax = vmax)
                    im1 = ax[1].imshow(pymap, extent = extent, aspect = 'auto', cmap = cmap, vmin = vmin, vmax = vmax)
                    for i, t in enumerate(['Px', 'Py']):
                        ax[i].invert_yscale()
                        ax[i].set_title(t)
                        ax[i].set_xlabel(D[0].get('Meta',{}).get('x_label', ''))
                        ax[i].set_ylabel(xlabel)
                if not type(pzmap) is list:
                    fig, ax = plt.subplots(figsize = figsize)
                    im2 = ax.imshow(pzmap, extent = extent, aspect = 'auto', cmap = cmap, vmin = vmin, vmax = vmax)
                    ax.invert_yscale()
                    ax.set_title('Pz')
                    ax.set_xlabel(D[0].get('Meta',{}).get('x_label', ''))
                    ax.set_ylabel(xlabel)
        #
        return results
    
    return results


# ==============================================================================================

def _polarization_xy(a_c1rp = None, a_c1rn = None, sf = None):
    """
    a_c1rp  : asymmetry from coil 1 and rotor +
    a_c1rn  : asymmetry from coil 1 and rotor -
    sf      : sherman function
    """
    if type(a_c1rp) is type(None) or type(a_c1rn) is type(None):
        return None, None
    if type(sf) is type(None):
        global SHERMAN
        sf = SHERMAN
    px = np.sqrt(2)/(2*sf)*(a_c1rp - a_c1rn)
    py = -np.sqrt(2)/(2*sf)*(a_c1rp + a_c1rn)
    return px, py

def _polarization_z(a_c2 = None, sf = None):
    """
    a_c1rp  : asymmetry from coil 1 and rotor +
    a_c1rn  : asymmetry from coil 1 and rotor -
    sf      : sherman function
    """
    if type(a_c2) is type(None):
        return None
    if type(sf) is type(None):
        global SHERMAN
        sf = SHERMAN
    pz = -a_c2/sf
    return pz







# ==============================================================================================
# ============================================================================================== F I T
# ==============================================================================================


# ----------------------------------- Gauss

def _gauss(x, *p):
    """
    p[0] * np.exp( -(x-p[1])**2 / (2*p[2]**2)) + p[3]*x + p[4]
    """
    A, mu, s, k, m = p
    return A * np.exp( -(x-mu)**2 / (2*s**2)) + k*x + m

def fitGauss(x, y, *p):
    """"""
    P = [None, None, None, None, None]
    if len(p) > 0: P[0] = p[0]  # A
    if len(p) > 1: P[1] = p[1]  # mu
    if len(p) > 2: P[2] = p[2]  # s
    if len(p) > 3: P[3] = p[3]  # k
    if len(p) > 4: P[4] = p[4]  # m
    #
    if type(P[3]) is type(None): # k
        P[3] = (y[-1] - y[0]) / (x[-1] - x[0])
    if type(P[4]) is type(None): # m
        P[4] = y[0] - P[3] * x[0]
    if type(P[0]) is type(None):
        i = abs(y - max(y)).argmin()
        P[0] = y[i] - (P[3] * x[i] + P[4])
    if type(P[1]) is type(None):
        i = abs(y - max(y)).argmin()
        P[1] = x[i]
    if type(P[2]) is type(None):
        P[2] = 0.01
    p = tuple(P)
    #print(", ".join(str(p) for p in P))
    #
    par, cov = curve_fit(_gauss, x, y, p0 = p)
    yfit = _gauss(x, *par)
    rdict = {'x': x,
             'y': y,
             'yfit': yfit,
             'par': par,
             'cov': cov,
             'profile': 'p0 * exp( -(x-p1)**2 / (2*p2**2)) + p3*x + p4',
             'method': _gauss}
    return rdict


def _gauss2(x, *p):
    """
    p0 * np.exp( -(x-p1)**2 / (2*p4**2)) + p2 * np.exp( -(x-p3)**2 / (2*p4**2)) + p5*x + p6
    """
    A1, mu1, A2, mu2, s, k, m = p
    return A1 * np.exp( -(x-mu1)**2 / (2*s**2)) + A2 * np.exp( -(x-mu2)**2 / (2*s**2)) + k*x + m

def fitGauss2(x, y, *p):
    """
    p0 * exp( -(x-p1)**2 / (2*p4**2)) + p2 * exp( -(x-p3)**2 / (2*p4**2)) + p5*x + p6
    """
    P = [None, None, None, None, None, None, None]
    if len(p) > 0: P[0] = p[0]  # A1
    if len(p) > 1: P[1] = p[1]  # mu1
    if len(p) > 2: P[2] = p[2]  # A2
    if len(p) > 3: P[3] = p[3]  # mu2
    if len(p) > 4: P[4] = p[4]  # s
    if len(p) > 5: P[5] = p[5]  # k
    if len(p) > 6: P[6] = p[6]  # m
    #
    if type(P[5]) is type(None):                # k
        P[5] = (y[-1] - y[0]) / (x[-1] - x[0])
    if type(P[6]) is type(None):                # m
        P[6] = y[0] - P[5] * x[0]
    if type(P[0]) is type(None):                # A1
        n = int(len(x)/2)
        i = abs(y[:n] - max(y[:n])).argmin()
        P[0] = y[i] - (P[5] * x[i] + P[5])
    if type(P[1]) is type(None):                # mu1
        i = abs(y[:n] - max(y[:n])).argmin()
        P[1] = x[i]
    if type(P[2]) is type(None):                # A2
        n = int(len(x)/2)
        i = abs(y[n:] - max(y[n:])).argmin()
        P[2] = y[n+i] - (P[5] * x[n+i] + P[6])
    if type(P[3]) is type(None):                # mu2
        i = abs(y[n:] - max(y[n:])).argmin()
        P[3] = x[n+i]
    if type(P[4]) is type(None):                # s
        P[4] = 0.01
    p = tuple(P)
    #print(", ".join(str(p) for p in P))
    #
    par, cov = curve_fit(_gauss2, x, y, p0 = p)
    yfit = _gauss2(x, *par)
    rdict = {'x': x,
             'y': y,
             'yfit': yfit,
             'par': par,
             'cov': cov,
             'profile': 'p0 * exp( -(x-p1)**2 / (2*p4**2)) + p2 * exp( -(x-p3)**2 / (2*p4**2)) + p5*x + p6',
             'method': _gauss2}
    return rdict



def PlotFitRes(D = {}, ax = None, **kwargs):
    """
    Plots data from a dict from Fit(), or fitGauss(), etc.

    Pyplot arguments: figsize, xlabel, ylabel, title, marker, s, linewidth, linestyle

    Additional arguments: only_fit, fit_density, marker_color, line_color, data_label, fit_label
    """
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must be a dict containing fit results.' + Fore.RESET); return
    if not('yfit' in D and 'par' in D and 'method' in D):
        print(Fore.RED + 'Argument D must be a dict containing fit results.' + Fore.RESET); return
    #
    recognized_kwargs = ['figsize', 'xlabel', 'ylabel', 'title', 'only_fit', 'marker', 's', 'fit_density',
                        'marker_color', 'line_color', 'linewidth', 'linestyle', 'data_label', 'fit_label']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)
    #
    figsize = kwargs.get('figsize', (5,2.5))
    xlabel = kwargs.get('xlabel', 'X')
    ylabel = kwargs.get('ylabel', 'Y')
    title = kwargs.get('title', 'fit')
    only_fit = kwargs.get('only_fit', False)
    marker = kwargs.get('marker', 'x')
    s = kwargs.get('s', 5)
    fit_density = kwargs.get('fit_density', 10)
    marker_color = kwargs.get('marker_color', 'tab:blue')
    line_color = kwargs.get('line_color', 'tab:orange')
    linewidth = kwargs.get('linewidth', 0.75)
    linestyle = kwargs.get('linestyle', '-') 
    data_label = kwargs.get('data_label', 'data')
    fit_label = kwargs.get('fit_label', 'fit')  
    #
    if not hasattr(ax, 'plot'):
        fig, ax = plt.subplots(figsize = figsize)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if not only_fit:
        ax.scatter(D['x'], D['y'], color = marker_color, marker = marker, s = s, label = data_label)
    xfit = np.linspace(D['x'][0], D['x'][-1], len(D['x']) * fit_density)
    yfit = D['method'](xfit, *D['par'])
    ax.plot(xfit, yfit, color = line_color, linewidth = linewidth, linestyle = linestyle, label = fit_label)
    ax.legend()
    return ax


def extendFit(fit = {}, x1 = None, x2 = None, n = None):
    """
    Extends the fit result from Fit().
    """
    if not type(fit) is dict:
        print(Fore.RED + 'Argument fit must be a dict from e.g. Fit() or similar method.' + Fore.RESET); return fit
    if not "method" in fit:
        print(Fore.RED + 'Argument fit must contain a reference to the fit profile used.' + Fore.RESET); return fit
    #
    xaxis = fit.get('x')
    y = fit.get('y')
    if type(xaxis) is type(None):
        print(Fore.RED + 'Argument fit has no key x.' + Fore.RESET); return fit
    if type(x1) is type(None): x1 = xaxis[0]
    if type(x2) is type(None): x2 = xaxis[-1]
    if type(n) is type(None): n = len(xaxis)
    #
    new_xaxis = np.linspace(x1, x2, n)
    new_yfit = fit['method'](new_xaxis, *tuple(fit['par']))
    new_y = np.copy(new_yfit) * np.NaN
    for i, X in enumerate(xaxis):
        indx = abs(X - new_xaxis).argmin()
        new_y[indx] = y[i]
    FIT = copy.deepcopy(fit)
    FIT.update({'x': new_xaxis})
    FIT.update({'y': new_y})
    FIT.update({'yfit': new_yfit})
    return FIT





def Fit(D = {}, profile = '', **kwargs):
    """
    Arguments:  profile     gauss, with argument p0 = [A,mu,s,k,m]
                            gauss2, with argument p0 = [A1,mu1,A2,mu2,s,k,m]
    """
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must be of type grumpy dict.' + Fore.RESET); return D
    profiles = ['gauss', 'gauss2']
    if profile == '':
        print('Pass argument profile as one of the following: {0}'.format(profiles))
        return D
    #
    recognized_kwargs = ['p0', 'y']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)
    #
    Int = D.get("int", np.array([]))
    dimension = len(np.shape(Int))
    #
    if dimension == 1:
        if profile == 'gauss':
            p0 = kwargs.get('p0', [None, None, None, None, None])
        if profile == 'gauss2':
            p0 = kwargs.get('p0', [None, None, None, None, None, None, None])
        if not 'asymmetry' in D:
            res = fitGauss(D['x'], D['int'], *tuple(p0))
        else:
            intkey = kwargs.get('y', '(none)')
            keys = ['int', 'int_on', 'asymmetry']
            if not intkey in keys:
                print(Fore.RED + 'For spin data, pass argument y as one of the following: {0}'.format(keys))
                return {}
            res = fitGauss(D['x'], D[intkey], *tuple(p0))
        return res
        
        # add more profiles here
    
    else: # continue with 2d fits here
        print(Fore.LIGHTRED_EX + "This method is currently only working for 1d data types." + Fore.RESET)


# ==============================================================================================
# ============================================================================================== S M O O T H
# ==============================================================================================



# ----------------------------------- Remove spikes

# This guy is (e.g.) used when loading spiky spind edc:s.
# Args to pass to Load are filter
def removeSpikes(Y,threshold=2):
    """
    Curtesy of Craig.
    Removes spikes from 1D or 2D intensity arrays.
    """
    if len(np.shape(Y)) > 2:
        print(Fore.RED + "This works for 1d and 2d data."+ Fore.RESET)
        return Y

    if len(np.shape(Y)) == 1:
        y_out = np.copy(Y)
        for i, y in enumerate(Y):
            if i > 0 and i < len(Y)-1:
                if y > (threshold*(Y[i-1] + Y[i+1])/2):
                    y_out[i]=(Y[i-1] + Y[i+1])/2
        return y_out
    
    if len(np.shape(Y)) == 2:
        YY = np.copy(Y)
        # the interior
        for iy in range(1, np.shape(Y)[0]-1):
            for ix in range(1, np.shape(Y)[1]-1):
                v = (YY[iy][ix-1] + YY[iy-1][ix] + YY[iy+1][ix] + YY[iy][ix+1])/4
                if YY[iy][ix] > threshold * v:
                    YY[iy][ix] = v
        # upper and lower edges
        for ix in range(1, np.shape(Y)[1]-1):
            v = (YY[0][ix-1] + YY[0][ix+1] + YY[1][ix])/3
            if YY[0][ix] > threshold * v:
                YY[0][ix] = v
            v = (YY[-1][ix-1] + YY[-1][ix+1] + YY[-2][ix])/3
            if YY[-1][ix] > threshold * v:
                YY[-1][ix] = v
        # left and right edges
        for iy in range(1, np.shape(Y)[0]-1):
            v = (YY[iy-1][0] + YY[iy+1][0] + YY[iy][1])/3
            if YY[iy][0] > threshold * v: 
                YY[iy][0] = v
            v = (YY[iy-1][-1] + YY[iy+1][-1] + YY[iy][-2])/3
            if YY[iy][-1] > threshold * v: 
                YY[iy][-1] = v
        # corners
        v = (YY[0][1] + YY[1][1] + YY[1][0])/3
        if YY[0][0] > threshold * v: YY[0][0] = v
        v = (YY[0][-2] + YY[1][-2] + YY[1][-1])/3
        if YY[0][-1] > threshold * v: YY[0][-1] = v
        v = (YY[-1][1] + YY[-2][1] + YY[-2][0])/3
        if YY[-1][0] > threshold * v: YY[-1][0] = v
        v = (YY[-1][-2] + YY[-2][-2] + YY[-2][-1])/3
        if YY[-1][-1] > threshold * v: YY[-1][-1] = v
        #
        return YY


def removeSpikes2(Y,threshold=2):
    """
    Removes spikes from 1D arrays.
    """
    if len(np.shape(Y)) > 1:
        print(Fore.RED + "This works for 1d data."+ Fore.RESET)
        return Y
    #
    for i, y in enumerate(Y[:-4]):
        ymin = min(Y[i:i+4])
        if Y[i] > threshold * ymin:
            Y[i] = (Y[i+1]+Y[i+2]+Y[i+3])/3
        if Y[i+1] > threshold * ymin:
            Y[i+1] = (Y[i]+Y[i+2]+Y[i+3])/3
        if Y[i+2] > threshold * ymin:
            Y[i+2] = (Y[i]+Y[i+1]+Y[i+3])/3
        if Y[i+3] > threshold * ymin:
            Y[i+3] = (Y[i]+Y[i+1]+Y[i+2])/3
    return Y
            
    
    



# median_filter from scipy.ndimage
def medianFilter(Y, size = 3, mode = 'reflect'):
    """
    Works for data of any dimention.
    
    Uses median_filter from scipy.ndimage
        scipy.ndimage.median_filter(input, size=None, footprint=None, 
                                    output=None, mode='reflect', cval=0.0, origin=0)
    kwargs:
        size, default 3
        mode, default 'reflect'. 
            can be 'reflect', 'constant', 'nearest', 'mirror', or 'wrap' where
                'reflect'   (d c b a | a b c d | d c b a)
                'constant'  (k k k k | a b c d | k k k k)
                'nearest'   (a a a a | a b c d | d d d d)
                'mirror'    (d c b | a b c d | c b a)
                'wrap'      (a b c d | a b c d | a b c d)

    """
    return median_filter(Y, size=size, mode=mode)
       



# a simple method for finding FWHM
def fwhm(x = [], y = [], x0 = None):
    """
    Pass x and y as equal-length arrays (or lists) and the method returns a dict with the FWHM, the
    intensity at HM, and the positions for HM (as well the indices). This is a crude method and works
    best for good statistics. Argument x0 can be used to help out when the simple algorithm screws up.
    For x0 = None the algorithm looks for y.max() and finds x1 and x2 at y.max()/2 on either side.
    """
    try: x, y = np.array(x), np.array(y)
    except:
        print(Fore.RED + "Arguments 'x' and 'y' must be lists or arrays."); return None, None
    if not len(x) == len(y):
        print(Fore.RED + "Arguments 'x' and 'y' must be lists or arrays of equal size."); return None, None
    if len(x) == 0:
        print(Fore.RED + "Arguments 'x' and 'y' can not have length zero."); return None, None
    if x[0] > x[-1]:
        x, y = x[::-1], y[::-1]
    #
    if type(x0) is type(None):
        indx0 = abs(y - y.max()).argmin()
        hm = (y.max() - y.min())/2
    else:
        try: x0 = float(x0)
        except:
            print(Fore.RED + "Argument 'x0' must be a number (or None to be ignored)." + Fore.RESET); return None, None
        if x0 < x.min() or x0 > x.max():
            print(Fore.RED + f"Argument 'x0' must be in the range ({x.min()}, {x.max()})." + Fore.RESET); return None, None
        indx0 = abs(x - x0).argmin()
        hm = (y[indx0] - y.min())/2
    indx1 = abs(y[:indx0] - hm).argmin()
    indx2 = abs(y[indx0:] - hm).argmin() + indx0
    res = {'fwhm': x[indx2] - x[indx1],
           'hm': hm,
           'x1': x[indx1],
           'x2': x[indx2],
           'i1': indx1,
           'i2': indx2}
    return res





def Smooth(D = {}, method = '', **kwargs):
    """
    Smooths intensity.
    Methods available now: 'median_filter', 'remove_spikes'
    Arguments:  for median filter: size and mode
                for remove_spikes: threshold
    Arguments for spin data: spin_key as 'int', 'int_on', or 'asymmetry'
    """
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must be a grumpy-type dict.' + Fore.RESET); return D
    if D == {}:
        print(Fore.RED + 'Argument D must be a grumpy-type dict.' + Fore.RESET); return D
    #
    recognized_kwargs = ['size', 'mode', 'threshold', 'spin']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)
    #
    mf_size = kwargs.get('size', 3)
    mf_mode = kwargs.get('mode', 'reflect')
    rs_threshold = kwargs.get('threshold', 1.2)
    spin_key = kwargs.get('spin', '')

    #
    methods = ['median_filter', 'remove_spikes']
    if method == '':
        print("Note: argument method was not passed so using default method = 'median_filter'.")
        method = 'median_filter'
    if not method in methods:
        msg = Fore.RESET + ", ".join(m for m in methods)
        print(Fore.RED + "Argument method must be one of: {0}.".format(msg) + Fore.RESET)
        return D
    #
    INTKEY = ''
    spin_keys = ['int', 'int_on', 'asymmetry']
    if 'asymmetry' in D: 
        if spin_key == '':
            print("Note: argument spin was not passed so using default spin = 'asymmetry'.")
            spin_key = 'asymmetry'
        if not spin_key in spin_keys:
            msg = Fore.RESET + ", ".join(m for m in methods)
            print(Fore.RED + "Argument method must be one of: {0}.".format(msg) + Fore.RESET)
            return D
        INTKEY = spin_key
    else:
        INTKEY = 'int'
    #
    dimension = len( np.shape( D.get(INTKEY, np.array([]))))
    if dimension == 1 or dimension == 2:
        if method == 'median_filter':
            Y = medianFilter(D.get(INTKEY), size = mf_size, mode = mf_mode)
        elif method == 'remove_spikes':
            Y = removeSpikes(D.get(INTKEY), threshold = rs_threshold)
    else:
        print(Fore.RED + "This method is only accepting 1d and 2d data." + Fore.RESET)
        return D
    #
    DD = copy.deepcopy(D)
    DD.update({INTKEY: Y})
    return DD

        
    
# ==============================================================================================
# ==============================================================================================
# ==============================================================================================
# ==============================================================================================

class MakeArithmetic():
    """
    Make a Grumpy dict into an object which you can perform arithmetic operation on, plus
    a few other stuff, like plotting. The axes and intensity are:
        self.x      (3d, 2d, 1d)
        self.y      (3d, 2d)
        self.z      (3d)
        self.int    (3d, 2d, 1d)
        
    Operators: 
        +, -, *, /, <, >, <=, >=
        
    Methods: 
        plot()  Plots 1d, and 2d data
                Optional arguments:
                    ax          as matplotlib.axes._axes.Axes
                    figsize     as tuple,eg (3,3) 
                    aspect      as "auto" or "equal"
                    vmin, vmax  as numbers or None
                    cbar        as boolean, adds a colorbar
                    cmap        as str or None, colormap e.g "bone_r", "bwr", "hot", "rainbow", ...
                    
        
    """
    def __init__(self, d = {}):
        """"""
        if not type(d) is dict:
            print(Fore.RED + "Argument d must be a Grumpy dict." + Fore.RESET)
            return
        #
        dm = 0
        axis = d.get("x", np.zeros([0]))
        if len(axis) > 0:
            self.x = axis
            dm += 1
        axis = d.get("y", np.zeros([0]))
        if len(axis) > 0:
            self.y = axis
            dm += 1
        axis = d.get("z", np.zeros([0]))
        if len(axis) > 0:
            self.z = axis
            dm += 1
        if dm == 0:
            print(Fore.RED + "Did not find any axes in the dict." + Fore.RESET)
        intensity = d.get("int", None)
        if type(intensity) is type(None):
            print(Fore.RED + "Did not find an intensity array in the dict." + Fore.RESET)
            return
        self.int = intensity
        if dm == 0: return
        #
        idim = len(np.shape(intensity))
        if not dm == idim:
            print(Fore.RED + f"The number of axes ({dm}) does not match the dimension of the intensity array ({idim})." + Fore.RESET)
            return
        #
        if dm == 1:
            if not len(self.x) == np.shape(intensity)[0]:
                print(Fore.RED + f"The length of the axis ({len(self.x)}) does not match the length of the intensity array ({np.shape(intensity)[0]})." + Fore.RESET)
                return
        elif dm == 2:
            if not len(self.x) == np.shape(intensity)[1] or not len(self.y) == np.shape(intensity)[0]:
                print(Fore.RED + f"The lengths of the axes ({len(self.y)}, {len(self.x)}) does not match the size of the intensity array ({np.shape(intensity)})." + Fore.RESET)
                return
        elif dm == 3:
            if not len(self.x) == np.shape(intensity)[2] or not len(self.y) == np.shape(intensity)[1] or not len(self.z) == np.shape(intensity)[0]:
                print(Fore.RED + f"The lengths of the axes ({len(self.z)},{len(self.y)}, {len(self.x)}) does not match the size of the intensity array ({np.shape(intensity)})." + Fore.RESET)
                return
    
    # --------------- operators -----------------
    
    def __add__(self, o):
        ret = copy.deepcopy(self)
        if type(o) is type(self):
            if self._check_other_object(o):
                ret.int = ret.int + o.int
        elif type(o) is int or type(o) is float:
            ret.int = ret.int + o
        else:
            print(Fore.RED + "The + operator only acts on same type, float, and int." + Fore.RESET)
        return ret
    
    def __sub__(self, o):
        ret = copy.deepcopy(self)
        if type(o) is type(self):
            if self._check_other_object(o):
                ret.int = ret.int - o.int
        elif type(o) is int or type(o) is float:
            ret.int = ret.int - o
        else:
            print(Fore.RED + "The - operator only acts on same type, float, and int." + Fore.RESET)
        return ret
    
    def __mul__(self, o):
        ret = copy.deepcopy(self)
        if type(o) is type(self):
            if self._check_other_object(o):
                ret.int = ret.int * o.int
        elif type(o) is int or type(o) is float:
            ret.int = ret.int * o
        else:
            print(Fore.RED + "The * operator only acts on same type (element-wise), float, and int." + Fore.RESET)
        return ret
    
    def __truediv__(self, o):
        ret = copy.deepcopy(self)
        if type(o) is type(self):
            if self._check_other_object(o):
                ret.int = ret.int / o.int
        elif type(o) is int or type(o) is float:
            ret.int = ret.int / o
        else:
            print(Fore.RED + "The / operator only acts on same type (element-wise), float, and int." + Fore.RESET)
        return ret
    
    def __lt__(self, o): # less than
        ret = False
        if type(o) is type(self):
            if self._check_other_object(o):
                ret = self.int.sum() < o.int.sum()
        elif type(o) is int or type(o) is float:
            ret = self.int.sum() < o
        else:
            print(Fore.RED + "The < operator only acts on same type, float, and int." + Fore.RESET)
        return ret
    
    def __gt__(self, o): # greater than
        ret = False
        if type(o) is type(self):
            if self._check_other_object(o):
                ret = self.int.sum() > o.int.sum()
        elif type(o) is int or type(o) is float:
            ret = self.int.sum() > o
        else:
            print(Fore.RED + "The > operator only acts on same type, float, and int." + Fore.RESET)
        return ret
    
    def __le__(self, o): # less or equal than
        ret = False
        if type(o) is type(self):
            if self._check_other_object(o):
                ret = self.int.sum() <= o.int.sum()
        elif type(o) is int or type(o) is float:
            ret = self.int.sum() <= o
        else:
            print(Fore.RED + "The <= operator only acts on same type, float, and int." + Fore.RESET)
        return ret
    
    def __ge__(self, o): # less or greater than
        ret = False
        if type(o) is type(self):
            if self._check_other_object(o):
                ret = self.int.sum() >= o.int.sum()
        elif type(o) is int or type(o) is float:
            ret = self.int.sum() >= o
        else:
            print(Fore.RED + "The >= operator only acts on same type, float, and int." + Fore.RESET)
        return ret
    
    def _check_other_object(self, o):
        try:
            oshape = np.shape(o.int)
        except:
            print(Fore.RED + "The data in the 2nd object is of unknown type." + Fore.RESET)
            return False
        if not oshape == np.shape(self.int):
            print(Fore.RED + "The data in the objects are not of the same size and/or dimentsion." + Fore.RESET)
            return False
        return True
    
    # --------------- plot -----------------
    
    def plot(self, ax = None, figsize = (3,3), aspect = "equal", vmin = None, vmax = None, cbar = False, cmap = None):
        try:
            idim = len(np.shape(self.int))
        except:
            print(Fore.RED + "There is no intensity data." + Fore.RESET)
            return ax
        if idim == 0:
            print(Fore.RED + "There is no intensity data." + Fore.RESET)
            return ax
        if not idim in [1, 2]:
            print(Fore.RED + "I'm not plotting data of higher dimension than 2." + Fore.RESET)
            return ax
        #
        if type(ax) is type(None):
            fig, ax = plt.subplots(figsize = figsize)
        #
        if idim == 1:
            ax.plot(self.x, self.int, linewidth = 0.75, color = "k")
            ax.set_xlabel("x"); ax.set_ylabel("intensity")
        elif idim == 2:
            extent = [self.y[0], self.y[-1], self.x[-1], self.x[0]]
            ims = ax.imshow(np.transpose(self.int), extent = extent, aspect = aspect, vmin = vmin, vmax = vmax, cmap = cmap)
            ax.invert_yaxis()
            if cbar:
                _ = plt.colorbar(ims, ax = ax)
        #
        return ax
    
    
            
            
    
        
        
            
            
        

            
            
            




# ==============================================================================================
# ==============================================================================================
# ==============================================================================================
# ==============================================================================================

def Dict(D = {}, **kwargs):
    """
    Just a method to show what is in a dict.
    """
    if not type(D) is dict: return
    if len(list(D.keys())) == 0: return
    #
    for key in list(D.keys()):
        item = D[key]
        if type(item) is dict:
            print(Fore.GREEN + f"{key}" + "{" + Fore.RESET)
            Dict(item)
            print(Fore.GREEN + '}' + Fore.RESET)
        if type(item) is str:
            print("{0:<22}{1}".format(key, item))
        if type(item) is list:
            print('{0:<22}list, {1}'.format(key, np.shape(item)))
        if type(item) is np.ndarray:
            print('{0:<22}array, {1}'.format(key, np.shape(item)))
        if type(item) is int or type(item) is float:
            print("{0:<22}{1}".format(key, item))




# ==============================================================================================
# ==============================================================================================
# ==============================================================================================
# ==============================================================================================

def _kwarg_checker(key_list = [], **kwargs):
    """
    Helper
    Reports non-existing kwargs.
    """
    if len(key_list) == 0 or len(kwargs) == 0: return
    if kwargs.get('kw', False):
        print(Fore.BLUE + 'Recognized **arguments: {0}'.format(", ".join(k for k in key_list)) + Fore.RESET)
    c = 0
    for kw in kwargs.keys():
        if not kw in key_list and not kw == 'kw':
            print(Fore.RED + "Argument {0} is not recognized.".format(kw) + Fore.RESET)
            c += 1
    if c > 0:
        print('Recognized **arguments: {0}'.format(", ".join(k for k in key_list)))





# ==============================================================================================
# ==============================================================================================
# ==============================================================================================
# ==============================================================================================

def SortP(d = {}, parameter_list = []):
    """Not ready. Pet project."""
    print('\nThis method is not ready.\n')
    if not type(d) is dict:
        print('Error: Argument d must be of type dict.'); return d
    RETDICT = {}
    PARS0 = d.get('Parameters', [])
    PARVS0 = d.get('Parameter_values', [])
    DATA = d.get('Data', np.array([]))
    if d.get('Experiment', {}) == {}:
        print("Error: Measurement info missing is in this dict."); return d
    if len(PARS0) == 0:
        print("Note: There are no parameters in this dict. Exit."); return d
    
    if not type(parameter_list) is list:
        print('Error: Argument parameter_list must be of type list. Yes, really.'); return d
    if len(parameter_list) == 0:
        print("Note: You have not specified any parameters.")
        print("      Available parameters are:")
        for i, p in enumerate(PARS0):
            print("      {0}: {1}".format(i, p))
        return d
    
    # --- check if the passed parameters exit
    plisterror = False
    for p in parameter_list:
        if not p in PARS0:
            print('Error: parameter {0} can not be found.'.format(p)); plisterror = True
    if plisterror:
        print('Exit'); return d
    
    # --- get parameters and values
    PARS = []; PARSV = []
    for p in parameter_list:
        indx = np.where(np.array(PARS0) == p)[0][0]
        PARS.append(PARS0[indx])
        PARSV.append(np.array(PARVS0[indx]))
    
    RETDICT.update({'File': d.get('File', {})})
    RETDICT.update({'Experiment': d.get('Experiment', {})})
    RETDICT.update({'Parameters': PARS})
    RETDICT.update({'Parameter_values': PARSV})

    PARVSu = []
    needed_num_cols = 1
    for pv in PARSV:
        tmp = np.unique(np.array(pv))
        PARVSu.append(tmp)
        needed_num_cols *= len(tmp)
    available_num_cols = len(d['Data'])
    
    if not needed_num_cols == available_num_cols:
        print('\nWarning: It might not be possible to sort the data based on the requested parameters.')
        print('         Required:')
        for i, pn in enumerate(PARS):
            print('         {0}: {1} {2} values'.format(i, pn, len(PARVSu[i])))
        print('         Total:     {0}'.format(needed_num_cols))
        print('         Available: {0}'.format(available_num_cols))
        print('         Then again... it might.')
    

    return RETDICT

# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================


def AppendFEmaps(data1 = {}, data2 = {}, shup = False):
    """
    This is a crude method of appending two Fermi maps. The assumption is that they have the same
    energy axes and analyzer slit angle axes, and that the shift-x axes are split in two regions.
    The dict in first argument should have lower shift-x values that that in dict2.
    """
    print("This method is under development.")
    if not (MeasurementType(data1, shup = True) == "fermi_map" and MeasurementType(data2, shup = True) == "fermi_map"):
        print("The data must be fermi maps"); return {}
    #
    try:
        x1, y1, z1, int1 = data1['x'], data1['y'], data1['z'], data1['int']
    except:
        print("could not get x, y, z, and int from data1"); return {}
    try:
        x2, y2, z2, int2 = data2['x'], data2['y'], data2['z'], data2['int']
    except:
        print("could not get x, y, z, and int from data2"); return {}
    if not(len(x1) == len(x2) and x1[0] == x2[0] and x1[-1] and x2[-1]):
        print("the x-axes (energy axes) are not the same."); return {}
    if not(len(y1) == len(y2) and y1[0] == y2[0] and y1[-1] and y2[-1]):
        print("the y-axes (analyzer slit angle axes) are not the same."); return {}
    #
    if z1[0] > z2[0]:
        if not shup:
            print(Fore.MAGENTA + "Note: is seems like the data in arg. data1 is starting at higher shift-x values")
            print("than the data in data2. This is not corrected for. Swap the argument unless you,")
            print("for some reason, want it to be like this." + Fore.RESET)
    if z1[-1] >= z2[0]:
        if not shup: 
            print(Fore.MAGENTA + "Note: is seems like there is an overlap of shift-x axis values. Use the")
            print("SubArray() method to remedy this (if you want)."+ Fore.RESET)
    #
    Z, int = [], []
    for i, z in enumerate(z1):
        Z.append(z)
        int.append(int1[i])
    for i, z in enumerate(z2):
        Z.append(z)
        int.append(int2[i])
    Z, int = np.array(Z), np.array(int)
    #
    result = {'Type': 'fermi_map', 'MeasurementType': 'fermi_map', 'Meta': data1.get("Meta", {})}
    result.update({'x': x1, 'y': y1, 'z': Z, 'int': int})
    return result






# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================

def ExportSpinEDC(data = {}, exclude = []):
    if not type(data) is dict:
        print(Fore.RED + "Argument data must be a grumpy dict" + Fore.RESET)
        return
    #
    if not type(exclude) is list:
        print(Fore.RESET + "Argument exclude must be a list of integers. Ignoring this argument." + Fore.RESET)
        exclude = []
    #
    if not data.get("Experiment", {}).get("Analyzer", "") == "PhoibosSpin":
        print(Fore.RED + "The dicts does not contain spin data." + Fore.RESET)
        return
    if not "NegativePolarity" in data.get("Parameters"):
        print(Fore.RED + "Expected 'NegativePolarity' as a parameter but could not find it." + Fore.RESET)
        return
    try:
        Parameter_values = data.get("Parameter_values", [])[0]
    except:
        print(Fore.RED + "Missing parameter values for 'NegativePolarity'." + Fore.RESET)
        return
    if len(Parameter_values) == 0:
        print(Fore.RED + "Missing parameter values for 'NegativePolarity'." + Fore.RESET)
        return
    #
    Data = data.get("Data", [])
    if not len(np.shape(Data)) == 3:
        print(Fore.RED + "The data in the dict is not formatted as I would expect." + Fore.RESET)
        return
    #
    fidn = data.get("Experiment",{}).get("Spectrum_ID")
    try: fidn = int(fidn)
    except:
        print(Fore.RED + "Could not find data spectrum ID" + Fore.RESET)
        return
    #
    file_name = f"id{fidn}.dat"
    f = open(file_name, "w")
    f.write(f"# File id: {fidn}\n")
    ret = data.get("Experiment",{}).get("Lens_Mode", "")
    if not ret == "": f.write(f"# Lens mode: {ret}\n")
    ret = data.get("Experiment",{}).get("Scan_Mode", "")
    if not ret == "": f.write(f"# Scan mode: {ret}\n")
    ret = data.get("Experiment",{}).get("Dwell_Time", "")
    if not ret == "": f.write(f"# Dwell time: {ret}\n")
    ret = data.get("Experiment",{}).get("Excitation_Energy", "")
    if not ret == "": f.write(f"# Excitation energy: {ret}\n")
    ret = data.get("Experiment",{}).get("Pass_Energy", "")
    if not ret == "": f.write(f"# Pass energy: {ret}\n")
    ret = data.get("Experiment",{}).get("Column_labels", "")
    if not ret == "": f.write(f"# Column labels: {ret}\n")
    #
    column_txt, excluded_txt = "# Columns: energy, OFF_avg, ON_avg", "# Excluded from the avgerages:"
    for i, par in enumerate(Parameter_values):
        column_txt = f"{column_txt}, {i}_{par.strip(' ')}"
        if i in exclude: excluded_txt = f"{excluded_txt} {i}_{par.strip(' ')}"
    f.write(f"{column_txt}\n")
    if len(exclude) > 0: f.write(f"{excluded_txt}\n")
    #
    energy = Data[0][0]
    avg_off, avg_on = np.zeros(len(energy)), np.zeros(len(energy))
    n_off, n_on = 0, 0
    for i, D in enumerate(Data):
        if not i in exclude:
            if "OFF" in Parameter_values[i]:
                avg_off += D[1]
                n_off += 1
            elif "ON" in Parameter_values[i]:
                avg_on += D[1]
                n_on += 1
    if n_off > 0: avg_off /= n_off
    if n_on >  0: avg_on  /= n_on
    #
    for ie, E in enumerate(energy):
        row = f"{E:7.3f}\t{avg_off[ie]:11.3f}\t{avg_on[ie]:11.3f}"
        for D in Data:
            row = f"{row}\t{D[1][ie]:11.3f}"
        f.write(f"{row}\n")
    #
    f.close()















# ============================================================================================================
# ============================================================================================================ Merge spin edcs
# ============================================================================================================
# ============================================================================================================

def MergeSpinEDC(data1 = {}, data2 = {}):
    """
    """
    if not type(data1) is dict or not type(data2) is dict:
        print(Fore.RED + f"Arguments data1 and data2 must be dicts." + Fore.RESET)
        return {}
    #
    if data2 == {}: return data1
    #
    if not (data1.get("Measurement_type") == "spin_edc" and data2.get("Measurement_type") == "spin_edc"):
        print(Fore.RED + f"Arguments data1 and data2 must be spin edc dicts." + Fore.RESET)
        return {}
    #
    shp1, shp2 = np.shape(data1.get("Data")), np.shape(data2.get("Data"))
    if not (shp1[1] == shp2[1] and shp1[2] == shp2[2]):
        print(Fore.RED + f"Arguments data1 and data2 must have the same energy axis." + Fore.RESET)
        return {}
    #
    if not(data1.get("Experiment", {}).get("Ep") == data2.get("Experiment", {}).get("Ep")):
        print(Fore.RED + f"Arguments data1 and data2 have different pass energies." + Fore.RESET)
        return {}
    #
    if not(data1.get("Data")[0][0][0] == data2.get("Data")[0][0][0] and data1.get("Data")[0][0][-1] == data2.get("Data")[0][0][-1]):
        print(Fore.RED + f"Arguments data1 and data2 must have the same energy axis." + Fore.RESET)
        return {}
    #
    File = deepcopy(data1.get("File", {}))
    File.update({"file_name": f'{data1.get("File", {}).get("file_name", "")}, {data2.get("File", {}).get("file_name", "")}'})
    File.update({"file": f'{data1.get("File", {}).get("file", "")}, {data2.get("File", {}).get("file", "")}'})
    Experiment = deepcopy(data1.get("Experiment", {}))
    Experiment.update({"Spectrum_ID": -1})
    Parameters = deepcopy(data1.get("Parameters", []))
    Parameter_values = []
    Parameter_values1, Parameter_values2 = data1.get("Parameter_values", []), data2.get("Parameter_values", [])
    for p in Parameter_values1[0]: Parameter_values.append(p)
    for p in Parameter_values2[0]: Parameter_values.append(p)
    Data = []
    for d in data1.get("Data", []): Data.append(d)
    for d in data2.get("Data", []): Data.append(d)
    Measurement_type = data1.get("Measurement_type", "")
    x = data1.get("x", np.array([]))
    
    intens1, intens2, n1, n2, intens1b, intens2b = np.zeros(len(x)), np.zeros(len(x)), 0, 0, [], []
    for i, p in enumerate(Parameter_values):
        
        if "OFF" in p:
            intens1 += Data[i][1]
            intens1b.append(Data[i][1])
            n1 += 1
        if "ON" in p:
            intens2 += Data[i][1]
            intens2b.append(Data[i][1])
            n2 += 1
    intens1, intens2 = np.array(intens1)/n1, np.array(intens2)/n2

    asymmetry = np.zeros([len(x)])
    for i in range(len(x)):
        denom = (intens1[i] + intens2[i])
        if not denom == 0:
            nom = (intens1[i] - intens2[i])
            asymmetry[i] = (nom/denom)
        else:
            asymmetry[i] = 0

    Meta, Type = data1.get("Meta", {}), data1.get("Type", "")

    NewDict = {"File": File, "Experiment":Experiment, "Parameters":Parameters, "Parameter_values":Parameter_values, 
                "Data": Data, "Measurement_type": Measurement_type, "x":x, "int": intens1, "int_on": intens2,
                "int_all": intens1b, "int_all_on": intens2b, "asymmetry": asymmetry, "Meta": Meta, "Type": Type}

    return NewDict




# ============================================================================================================ Save spin EDCs to file

def SaveSpinEDC2File(data = [], file_name = "spin_edc.dat"):
    """
    Saves spin edc's to two text files (for negative polarity on and off) as columns, the first column being the energy.
    Arguments:
        data        A spin edc dict or a list of spin edc dicts.
        file_name   The file name as a string. Eg. 'spin.dat' will generate two files: spin_off.dat and spin_on.dat.
    """
    if type(data) is dict: data = [data]
    if not type(data) is list:
        print(Fore.Red + "The argument data must be a list of (at least 1) spin edc dicts." + Fore.RESET); return
    if len(data) == 0:
        print(Fore.Red + "The argument data must be a list of (at least 1) spin edc dicts." + Fore.RESET); return
    for i, d in enumerate(data):
        if not type(d) is dict:
            print(Fore.Red + f"Item {i+1} in argument data is not a dict." + Fore.RESET); return
        if not d.get("Measurement_type") == "spin_edc":
            print(Fore.Red + f"Item {i+1} in argument data is not spin edc dict." + Fore.RESET); return
    if not type(file_name) is str:
        print(Fore.MAGENTA + "The argument file_name must be a str. Setting it to default 'spin_edc.dat'" + Fore.RESET)
        file_name = "spin_edc.dat"
    #
    energy_axes = []
    for d in data: energy_axes.append(d["x"])
    e1, e2, el = energy_axes[0][0], energy_axes[0][-1], len(energy_axes[0])
    for e in energy_axes:
        if not (e[0] == e1 and e[-1] == e2 and len(e) == el):
            print(Fore.Red + "Not all spin edc's have the same energy axis." + Fore.RESET); return
    #
    spectrum_id = []
    for d in data:
        spectrum_id.append(d.get("Experiment", {}).get("Spectrum_ID", "?"))
    #
    lens_modes, pass_energies = [], []
    for d in data:
        lens_modes.append(d.get("Experiment", {}).get("Lens_Mode", "?"))
        pass_energies.append(d.get("Experiment", {}).get("Ep", "?"))
    #
    file_name_off = file_name.split(".")[0] + "_off." + file_name.split(".")[1]
    file_name_on  = file_name.split(".")[0] + "_on."  + file_name.split(".")[1]
    foff, fon = open(file_name_off, "w"), open(file_name_on, "w")
    #
    if len(spectrum_id) == 1:
        if int(spectrum_id[0]) > 0:
            foff.write(f"# spectrum id: {int(spectrum_id[0])}\n")
            fon.write(f"# spectrum id: {int(spectrum_id[0])}\n")
        else:
            foff.write(f"# spectrum id: several scans combined into one.\n")
            fon.write(f"# spectrum id: several scans combined into one.\n")
    else:
        foff.write(f"# spectrum id: {spectrum_id}\n")
        fon.write(f"# spectrum id: {spectrum_id}\n") 
    #
    foff.write("# negative polarity: off\n")
    fon.write("# negative polarity: on\n")
    #
    if len(lens_modes) == 1:
        foff.write(f"# lens mode: {lens_modes[0]}\n")
        fon.write(f"# lens mode: {lens_modes[0]}\n")
    else:
        if all(element == lens_modes[0] for element in lens_modes):
            foff.write(f"# lens mode: {lens_modes[0]}\n")
            fon.write(f"# lens mode: {lens_modes[0]}\n")
        else:
            print(Fore.MAGENTA + "Note: Different lens_modes were used." + Fore.RESET)
            foff.write(f"# lens mode: {lens_modes}\n")
            fon.write(f"# lens mode: {lens_modes}\n")
    if len(pass_energies) == 1:
        foff.write(f"# pass energy: {pass_energies[0]}\n")
        fon.write(f"# pass energy: {pass_energies[0]}\n")
    else:
        if all(element == pass_energies[0] for element in pass_energies):
            foff.write(f"# pass energy: {pass_energies[0]}\n")
            fon.write(f"# pass energy: {pass_energies[0]}\n")
        else:
            print(Fore.MAGENTA + "Note: Different pass energies were used." + Fore.RESET)
            foff.write(f"# pass energy: {pass_energies}\n")
            fon.write(f"# pass energy: {pass_energies}\n")
    del lens_modes, pass_energies
    #
    ENERGY = data[0]["x"]
    CURVES_OFF, CURVES_ON = [], []
    for D in data:
        for curve in D["int_all"]: CURVES_OFF.append(curve)
        for curve in D["int_all_on"]: CURVES_ON.append(curve)

    CURVES_OFF, CURVES_ON = np.array(CURVES_OFF), np.array(CURVES_ON)
    for i in range(len(ENERGY)-1):
        rowoff, rowon = f"{ENERGY[i]:.2f}", f"{ENERGY[i]:.2f}"
        for j in range(len(CURVES_OFF)):
            rowoff = f"{rowoff}\t{CURVES_OFF[j][i]:.1f}"
        for j in range(len(CURVES_ON)):
            rowon = f"{rowon}\t{CURVES_ON[j][i]:.1f}"
        foff.write(f"{rowoff}\n")
        fon.write(f"{rowon}\n")
    foff.close()
    fon.close()
    print(Fore.BLUE + f"Data saved to {file_name_off} and {file_name_on}." + Fore.RESET)





# ============================================================================================================ 

def RemoveSpinEDC(data = {}, edc_number = -1):
    """
    Remove a specific curve from the data. Returns a new data dict.
    Arguments:
        data        Spin edc dict

    """
    if not type(data) is dict:
        print(Fore.RED + "Argument data must be a spin edc dict." + Fore.RESET); return data
    if not data.get("Measurement_type") == "spin_edc":
        print(Fore.RED + "Argument data must be a spin edc dict." + Fore.RESET); return data
    try: edc_number = int(edc_number)
    except:
        print(Fore.RED + "Argument edc_number must be an integer." + Fore.RESET); return data
    num_curves = len(data.get("Data", -1))
    if edc_number >= num_curves or edc_number < 0:
        print(Fore.RED + f"Argument edc_number must be an integer from 0 to {num_curves - 1}." + Fore.RESET); return data
    #
    Data, Parameter_values, int_all, int_all_on = [],  [], [], []
    intensity = np.zeros([len(data["x"])])
    intensity_on = np.zeros([len(data["x"])])
    ion, ioff = 0, 0
    for i, curve in enumerate(data["Data"]):
        if not i == edc_number:
            Data.append(curve)
            Parameter_values.append(data["Parameter_values"][0][i])
            if "OFF" in Parameter_values[-1]:
                int_all.append(curve[1])
                intensity += curve[1]
                ioff += 1
            else:
                int_all_on.append(curve[1])
                intensity_on += curve[1]
                ion += 1
    if ioff > 0: intensity = intensity/ioff
    if ion > 0: intensity_on = intensity_on/ion
    #
    asymmetry = np.zeros([len(intensity)])
    for i in range(len(intensity)):
        denom = (intensity[i] + intensity_on[i])
        if not denom == 0:
            nom = (intensity[i] - intensity_on[i])
            asymmetry[i] = (nom/denom)
        else:
            asymmetry[i] = 0
    #
    newDict = deepcopy(data)
    newDict.update({"Data": Data})
    newDict.update({"Parameter_values": [Parameter_values]})
    newDict.update({"int": intensity})
    newDict.update({"int_on": intensity_on})
    newDict.update({"int_all": int_all})
    newDict.update({"int_all_on": int_all_on})
    newDict.update({"asymmetry": asymmetry})
    return newDict





# ============================================================================================================

class MassageSpinEDC():
    """
    """
    def __init__(self, data = {}):
        print(Fore.MAGENTA + "Note: I haven't been able to figure out how to make the decrease and increase buttons")
        print("automatically update the graph. For now: update by sliding a slider forwards and backwards.\n" + Fore.RESET)
        #
        print(Fore.BLUE + "Select curve, position the energy at a spike. Click the decrease and increase buttons")
        print("to adjust. Use the result() method to get an updated data dict.\n" + Fore.RESET)
        try:
            mtype = data.get("Measurement_type", "")
        except:
            print(Fore.RED + 'Argument data must be a grumpy dict.' + Fore.RESET); return
        if not mtype == "spin_edc":
            print(Fore.RED + 'Argument data must be a grumpy dict containing spin edc data.' + Fore.RESET); return
        #
        self.data = data
        self.Data = data.get("Data")
        self.energy = self.data.get("x", np.array([]))
        dX = (self.energy[-1]-self.energy[0])/(len(self.energy)-1)

        self._SliderC = ipw.FloatSlider(min=1, max=len(self.Data), step = 1, description = 'Curve', value = 0, readout_format = ".0f")
        self._SliderX = ipw.FloatSlider(min=self.energy[0], max=self.energy[-1], step = dX, description = 'Energy', value = dX, readout_format = ".2f")
        self.interact_args = {"curve": self._SliderC, "energy": self._SliderX}
        self._box1 = ipw.VBox([self._SliderC, self._SliderX ])

        self._ButtonPlus = ipw.Button(descripton = "Increase")
        self._ButtonMinus = ipw.Button(descripton = "Decrease")
        self._box2 = ipw.VBox([self._ButtonPlus, self._ButtonMinus])
        self._ButtonPlus.on_click(self._ButtonPlusCicked)
        self._ButtonMinus.on_click(self._ButtonMinusCicked)

        self.markers_checkbox = ipw.Checkbox(value = False, description = "Markers")
        self.legend_checkbox = ipw.Checkbox(value = False, description = "Legend")
        self.interact_args.update({"markers": self.markers_checkbox, "legend": self.legend_checkbox})
        self._box3 = ipw.VBox([self.markers_checkbox, self.legend_checkbox])

        def plot(curve, energy, markers, legend):#, button_plus, button_minus):
            Y = self.Data[int(curve)-1][1]
            fig, ax = plt.subplots(figsize = (7,4))

            if "OFF" in data.get("Parameter_values")[0][int(curve)-1]:
                lcol = "blue"
                lbl = f"curve {curve:.0f}, off"
            else: 
                lcol = "red"
                lbl = f"curve {curve:.0f}, on"
            if self.markers_checkbox.value:
                ax.plot(self.energy, Y, linestyle = "--", color = lcol, marker = "x", markersize = 5, linewidth = 0.6, label = lbl)
            else:
                ax.plot(self.energy, Y, linestyle = "-", color = lcol, label = lbl)

            ax.plot(self.energy, Y, color = lcol)
            ax.axvline(x = energy, color = 'k', linestyle = "--", linewidth = 0.6)
            ix = abs(energy - self.energy).argmin()
            ax.axhline(y = Y[ix], color = 'k', linestyle = "--", linewidth = 0.6)

            if self.legend_checkbox.value: ax.legend()

        self._Interact = ipw.interactive_output(plot, self.interact_args)
        self._box_out = ipw.VBox([ipw.HBox([self._box1, self._box2, self._box3]), self._Interact])
        display(self._box_out)

    def _ButtonPlusCicked(self, b):
        c = self._SliderC.value
        x = self._SliderX.value
        ix = abs(x - self.Data[int(c)-1][0]).argmin()
        y = self.Data[int(c)-1][1][ix]
        self.Data[int(c)-1][1][ix] = 1.025*y
        return 1

    def _ButtonMinusCicked(self, b):
        c = self._SliderC.value
        x = self._SliderX.value
        ix = abs(x - self.Data[int(c)-1][0]).argmin()
        y = self.Data[int(c)-1][1][ix]
        self.Data[int(c)-1][1][ix] = 0.975*y
        return 1


    def result(self):
        result = deepcopy(self.data)
        Parameter_values = result.get("Parameter_values", [])
        intens1, intens2, n1, n2, intens1b, intens2b = np.zeros(len(self.energy)), np.zeros(len(self.energy)), 0, 0, [], []
        for i, p in enumerate(Parameter_values[0]):
            if "OFF" in p:
                intens1 += self.Data[i][1]
                intens1b.append(self.Data[i][1])
                n1 += 1
            if "ON" in p:
                intens2 += self.Data[i][1]
                intens2b.append(self.Data[i][1])
                n2 += 1
        intens1, intens2 = np.array(intens1)/n1, np.array(intens2)/n2
        asymmetry = np.zeros([len(self.energy)])
        for i in range(len(self.energy)):
            denom = (intens1[i] + intens2[i])
            if not denom == 0:
                nom = (intens1[i] - intens2[i])
                asymmetry[i] = (nom/denom)
            else:
                asymmetry[i] = 0
        global SHERMAN
        result.update({"Data": self.Data, "int": intens1, "int_on": intens2, 
            "int_all": intens1b, "int_all_on": intens2b, "asymmetry": asymmetry/SHERMAN})
        return result







# ============================================================================================================
# ============================================================================================================ S A V E to T E X T
# ============================================================================================================

def Save(data = {}, file_name = "data.dat", **kwargs):
    """
    """
    if not type(data) is dict:
        print(Fore.RED + "Argument data must be a dict." + Fore.RESET); return
    if not type(file_name) is str:
        print(Fore.RED + "Argument file_name must be a string." + Fore.RESET); return
    #
    Measurement_type = data.get("Measurement_type", "")
    
    # ----------------------- spin edc
    if Measurement_type == "spin_edc":
        SaveSpinEDC2File(data = data, file_name = file_name)
        return

    # ----------------------- arpes
    elif Measurement_type == "ARPES":
        SaveARPES2File(data = data, file_name = file_name)
    
    # ----------------------- fermi map
    elif Measurement_type == "fermi_map":
        SaveFermiMap2File(data = data, file_name = file_name)
    
    # ----------------------- XPS
    elif Measurement_type == "XPS":
        SaveXPS2File(data = data, file_name = file_name)
    
    # ----------------------- scattering spectra
    elif Measurement_type == "target_scattering_spectrum":
        SaveTargetScatteringSpectrum2File(data = data, file_name = file_name)
    
    # ----------------------- else
    else:
        print(Fore.MAGENTA + "It appears that I have not implemented a method for saving this type of data to file yet. Stay tuned for updates." + Fore.RESET)


def SaveTargetScatteringSpectrum2File(data = {}, file_name = "data.dat"):
    if not type(data) is dict:
        print(Fore.RED + "Argument data must be a dict." + Fore.RESET); return
    if not type(file_name) is str:
        print(Fore.RED + "Argument file_name must be a string." + Fore.RESET); return
    if not data.get("Measurement_type", "") == "XPS":
        print(Fore.RED + "Argument data does not contain XPS data." + Fore.RESET); return
    #
    f = open(file_name, "w")
    f.write("# data: target scattering spectrum\n")
    f.write(f"# label x: {data.get('Meta', {}).get('x_label', '?')}\n")
    f.write(f"# label y: {data.get('Meta', {}).get('int_label', '?')}\n")
    f.write(f"# spectrum id: {data.get('Experiment', {}).get('Spectrum_ID', '?')}\n")
    f.write(f"# lens mode: {data.get('Experiment', {}).get('Lens_Mode', '?')}\n")
    f.write(f"# scan mode: {data.get('Experiment', {}).get('Scan_Mode', '?')}\n")
    f.write(f"# dwell time: {data.get('Experiment', {}).get('Dwell_Time', '?')}\n")
    f.write(f"# photon energy: {data.get('Experiment', {}).get('Excitation_Energy', '?')}\n")
    f.write(f"# kinetic energy: {data.get('Experiment', {}).get('Kinetic_Energy', '?')}\n")
    f.write(f"# pass energy: {data.get('Experiment', {}).get('Pass_Energy', '?')}\n")
    for i, x in data.get("x"): f.write(f"{x:7.3f}\t{data.get('int')[i]:.5e}\n")
    f.close()




def SaveXPS2File(data = {}, file_name = "data.dat"):
    if not type(data) is dict:
        print(Fore.RED + "Argument data must be a dict." + Fore.RESET); return
    if not type(file_name) is str:
        print(Fore.RED + "Argument file_name must be a string." + Fore.RESET); return
    if not data.get("Measurement_type", "") == "XPS":
        print(Fore.RED + "Argument data does not contain XPS data." + Fore.RESET); return
    #
    f = open(file_name, "w")
    f.write("# data: XPS\n")
    f.write(f"# label x: {data.get('Meta', {}).get('x_label', '?')}\n")
    f.write(f"# label y: {data.get('Meta', {}).get('int_label', '?')}\n")
    f.write(f"# spectrum id: {data.get('Experiment', {}).get('Spectrum_ID', '?')}\n")
    f.write(f"# lens mode: {data.get('Experiment', {}).get('Lens_Mode', '?')}\n")
    f.write(f"# scan mode: {data.get('Experiment', {}).get('Scan_Mode', '?')}\n")
    f.write(f"# dwell time: {data.get('Experiment', {}).get('Dwell_Time', '?')}\n")
    f.write(f"# photon energy: {data.get('Experiment', {}).get('Excitation_Energy', '?')}\n")
    f.write(f"# kinetic energy: {data.get('Experiment', {}).get('Kinetic_Energy', '?')}\n")
    f.write(f"# pass energy: {data.get('Experiment', {}).get('Pass_Energy', '?')}\n")
    for i, x in data.get("x"): f.write(f"{x:7.3f}\t{data.get('int')[i]:.5e}\n")
    f.close()



def SaveFermiMap2File(data = {}, file_name = "data.dat"):
    if not type(data) is dict:
        print(Fore.RED + "Argument data must be a dict." + Fore.RESET); return
    if not type(file_name) is str:
        print(Fore.RED + "Argument file_name must be a string." + Fore.RESET); return
    if not data.get("Measurement_type", "") == "fermi_map":
        print(Fore.RED + "Argument data does not contain Fermi map data." + Fore.RESET); return
    #
    fn_e = file_name.split(".")[0] + "_energy." + file_name.split(".")[1]
    fn_a = file_name.split(".")[0] + "_angle." + file_name.split(".")[1]
    fn_x = file_name.split(".")[0] + "_shiftx." + file_name.split(".")[1]
    fn_i = file_name.split(".")[0] + "_intensity." + file_name.split(".")[1]
    f_e, f_a, f_x, f_i = open(fn_e, "w"), open(fn_a, "w"), open(fn_x, "w"), open(fn_i, "w")
    F = [f_e, f_a, f_x, f_i]
    #
    F[0].write("# data: Fermi map, energy axis\n")
    F[0].write(f"# label: {data.get('Meta', {}).get('x_label', '?')}\n")
    F[1].write("# data: Fermi map, angle axis\n")
    F[1].write(f"# label: {data.get('Meta', {}).get('y_label', '?')}\n")
    F[2].write("# data: Fermi map, ShiftX\n")
    F[2].write(f"# label: {data.get('Meta', {}).get('z_label', '?')}\n")
    F[2].write("# data: Fermi map, intensity\n")
    F[2].write(f"# label: {data.get('Meta', {}).get('int_label', '?')}\n")
    for f in F: f.write(f"# spectrum id: {data.get('Experiment', {}).get('Spectrum_ID', '?')}\n")
    for f in F: f.write(f"# lens mode: {data.get('Experiment', {}).get('Lens_Mode', '?')}\n")
    for f in F: f.write(f"# scan mode: {data.get('Experiment', {}).get('Scan_Mode', '?')}\n")
    for f in F: f.write(f"# dwell time: {data.get('Experiment', {}).get('Dwell_Time', '?')}\n")
    for f in F: f.write(f"# photon energy: {data.get('Experiment', {}).get('Excitation_Energy', '?')}\n")
    for f in F: f.write(f"# kinetic energy: {data.get('Experiment', {}).get('Kinetic_Energy', '?')}\n")
    for f in F: f.write(f"# pass energy: {data.get('Experiment', {}).get('Pass_Energy', '?')}\n")
    F[3].write("#\n")
    lx, ly, le = len(data.get("z")), len(data.get("y")), len(data.get("x"))
    F[3].write(f"# The data is stored in one column of length {lx*ly*le}. Reshape it into \n")
    F[3].write(f"# a 3d array of shape (ShiftX = {lx}, SlitY = {ly}, E = {le}).\n#\n")
    #
    for x in data.get("x"): F[0].write(f"{x:7.3f}\n")
    F[0].close()
    for y in data.get("y"): F[1].write(f"{y:7.3f}\n")
    F[1].close()
    for z in data.get("z"): F[2].write(f"{z:7.3f}\n")
    F[2].close()
    #
    for i in range(len(data.get("z"))):
        F[3].write(f"# slice {i}, ShiftX = {data.get("z")[i]}\n")
        np.savetxt(F[3], data.get("int")[i,:,:].flatten(), fmt='%-.5e')
    F[3].close()






def SaveARPES2File(data = {}, file_name = "data.dat"):
    if not type(data) is dict:
        print(Fore.RED + "Argument data must be a dict." + Fore.RESET); return
    if not type(file_name) is str:
        print(Fore.RED + "Argument file_name must be a string." + Fore.RESET); return
    if not data.get("Measurement_type", "") == "ARPES":
        print(Fore.RED + "Argument data does not contain ARPES data." + Fore.RESET); return
    #
    fn_e = file_name.split(".")[0] + "_energy." + file_name.split(".")[1]
    fn_a = file_name.split(".")[0] + "_angle." + file_name.split(".")[1]
    fn_i = file_name.split(".")[0] + "_intensity." + file_name.split(".")[1]
    f_e, f_a, f_i = open(fn_e, "w"), open(fn_a, "w"), open(fn_i, "w")
    F = [f_e, f_a, f_i]
    #
    F[0].write("# data: ARPES, energy axis\n")
    F[0].write(f"# label: {data.get('Meta', {}).get('x_label', '?')}\n")
    F[1].write("# data: ARPES, angle axis\n")
    F[1].write(f"# label: {data.get('Meta', {}).get('y_label', '?')}\n")
    F[2].write("# data: ARPES, intensity\n")
    F[2].write(f"# label: {data.get('Meta', {}).get('int_label', '?')}\n")
    for f in F: f.write(f"# spectrum id: {data.get('Experiment', {}).get('Spectrum_ID', '?')}\n")
    for f in F: f.write(f"# lens mode: {data.get('Experiment', {}).get('Lens_Mode', '?')}\n")
    for f in F: f.write(f"# scan mode: {data.get('Experiment', {}).get('Scan_Mode', '?')}\n")
    for f in F: f.write(f"# dwell time: {data.get('Experiment', {}).get('Dwell_Time', '?')}\n")
    for f in F: f.write(f"# photon energy: {data.get('Experiment', {}).get('Excitation_Energy', '?')}\n")
    for f in F: f.write(f"# kinetic energy: {data.get('Experiment', {}).get('Kinetic_Energy', '?')}\n")
    for f in F: f.write(f"# pass energy: {data.get('Experiment', {}).get('Pass_Energy', '?')}\n")
    #
    for x in data.get("x"): F[0].write(f"{x:7.3f}\n")
    F[0].close()
    for y in data.get("y"): F[1].write(f"{y:7.3f}\n")
    F[1].close()
    #
    for I in data.get("int"):
        row = f"{I[0]:.5e}"
        for index, i in enumerate(I):
            if index> 0: row = f"{row}\t{i:.5e}"
        F[2].write(f"{row}\n")
    F[2].close()









