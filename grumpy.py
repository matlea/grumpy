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
    QuickSpinMDC()
    Polarization()      Not checked after re-write.

    Explore()           For 2d and 3d data. Not quite ready but it works. Need to add option to rotate the images.
    ShiftAxis()
    PlotFermiMap()      (Is called by Plot() when needed)
    Fit()               Work in progress. Ready for 1d data with profiles gauss and gauss2.
    Smooth()            Work in progress. Ready for 1d and 2d data.
    PlotFitRes()        Plots fit results (also called by Plot() if needed).
    info()
    Info()

    AppendFEmaps()      An unsophiticated method to append two Fermi maps.

    xps()               Not a method but a class, containing binding energies for elements.

    SortP()             (Loader for general data not covered by the types listed above. A work in progress. Slow work, slow progress.)

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

__version__ = "23.09.27"
__author__  = "Mats Leandersson"


import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit

from colorama import Fore
from sys import getsizeof

try: 
    import ipywidgets as ipw
    from IPython.display import display
except: 
    print(Fore.RED + '(\nCould not import the ipywidget module and/or display from IPython.display.') 
    print('Interactive plots will not work.\n\n' + Fore.BLACK)

# ==============================================================================================
print(Fore.LIGHTWHITE_EX + f"--- grumpy, version {__version__} (SAL X = ShiftX, SALY = -ShiftY, Asymmetry: negative polarity off - on.)" + Fore.BLACK)
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





# ============================================================================================================
# ============================================================================================================
# ============================================================================================================
# ============================================================================================================


def MeasurementType(data, shup = False, **kwargs):
    if not type(data) is dict:
        print(Fore.RED + 'Argument data must be a grumpy dict.' + Fore.BLACK); return
    
    Experiment = data.get('Experiment', {})
    
    Measurement_type = ''
    #global Measurement_types
    
    ANALYZER = Experiment.get('Analyzer', '')
    NON_ENERGY_ORDINATE = np.array(data.get('Non_Energy_Ordinate', np.array([])))
    PARAMETERS = data.get('Parameters', [])

    # ARPES or XPS
    if len(PARAMETERS) == 0 and ANALYZER == 'PhoibosCCD' and len(NON_ENERGY_ORDINATE) > 0:
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
    elif len(PARAMETERS) == 1 and ("SAL X [deg]" in PARAMETERS or "ShiftX [a.u.]" in PARAMETERS) and ANALYZER == 'PhoibosCCD':
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



def Load(file_name = '', path = '', shup = False, **kwargs):
    """
    Load exported data from Prodigy.
    Identifies common measurement setups.
    """
    recognized_kwargs = ['remove_spikes', 'filter_outliers', 'threshold']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    DC = {}             # The main dict
    File = {}           # Sub dict, add to dict DC
    Experiment = {}     # Sub dict, add to dict DC
    Parameters = []     # Array, add to dict Experiment

    try: nshup = not shup
    except: nshup = False

    if not type(file_name) is str or not type(path) is str:
        print(Fore.RED + "Args file_name and path must be strings." + Fore.BLACK)
        return {}
    if file_name == '':
        print(Fore.RED + "Arg file_name can not be an empty string." + Fore.BLACK)
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
        print(Fore.RED + 'Could not find/open the file.' + Fore.BLACK)
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
                        Data.append([float(values[0]), float(values[-1])])
                    except:
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
        if abs(angle_y.min()) < 1 and abs(angle_y.max()) < 1: 
            DC.update({'Measurement_type': Measurement_type})
            DC.update({'x': energy})
            Meta.update({'x_label': DC.get('Experiment').get('Energy_Axis', '')})
            Meta.update({'x_type': defax_energy})
            DC.update({'int': data[:,:].sum(axis = 0)})
            Meta.update({'intlabel': DC.get('Experiment').get('Count_Rate', '')})
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
                print(Fore.LIGHTRED_EX + "(Note: in the future use arg. remove_spikes instead of filter_outliers.)"+ Fore.BLACK)
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
                print(Fore.LIGHTRED_EX + "(Note: in the future use arg. remove_spikes instead of filter_outliers.)"+ Fore.BLACK)
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
                print(Fore.LIGHTRED_EX + "(Note: in the future use arg. remove_spikes instead of filter_outliers.)"+ Fore.BLACK)

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
    RED = Fore.RED; BLU = Fore.BLUE; BLK = Fore.BLACK
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
    BLK = Fore.BLACK
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
    print('int:                {0}'.format(shp) + Fore.BLACK)
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
        print(Fore.RED + 'Argument D must be a grumpy dict.' + Fore.BLACK); return D
    #
    av_axes = _getAxes(D)
    if len(av_axes) == 0:
        print(Fore.RED + 'Could not find any axes in the dict.' + Fore.BLACK); return D
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
        print(Fore.RED + "There is no intensity data in the dict." + Fore.BLACK); return D
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
            print(Fore.RED + 'Deflector spin not implemented yet.' + Fore.BLACK); return D
            
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
        print(Fore.RED + 'Argument D must be a grumpy dict.'+ Fore.BLACK); return D
    #
    av_axes = _getAxes(D)
    if len(av_axes) == 0:
        print(Fore.RED + 'Could not find any axes in the dict.' + Fore.BLACK); return D
    #
    if D.get('Type').startswith('spin'):
        print(Fore.RED + 'Compact does not work with spin data yet.' + Fore.BLACK); return D
    #
    axis = axis.lower()
    if not axis in list(np.transpose(av_axes)[0]):
        print(Fore.RED + "Axis {0} is not found in the dict.".format(axis) + Fore.BLACK); return D
    #
    Int = D.get('int')
    if type(Int) is None:
        print(Fore.RED + "There is no intensity data in the dict." + Fore.BLACK); return D
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
        print(Fore.RED + 'Argument D must be a grumpy dict.' + Fore.BLACK); return D
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
        print(Fore.RED + 'Can not find any axes in argument D.' + Fore.BLACK); return D
    axes = list(np.transpose(axes_info))[0]
    #
    if axis == '':
        print(Fore.RED + 'Argument axis must be one of the following: {0}'.format(axes) + Fore.BLACK); return D
    #
    if not axis in axes:
        print(Fore.RED + 'Argument axis must be one of the following: {0}'.format(axes) + Fore.BLACK); return D
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
        print(Fore.RED + 'Argument D must be a grumpy dict.' + Fore.BLACK); return D
    #
    axes_info = _getAxes(D = D)
    axes = list(np.transpose(axes_info))[0]
    if not axis in axes:
        print(Fore.RED + "axis '{0}' is not found in this grumpy dict.".format(axis) + Fore.BLACK); return D
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
        print(Fore.RED + 'Argument D must be a grumpy dict.'+ Fore.BLACK); return ax

    
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
        print(Fore.RED + 'The data is not recognized by this method. Plot it yourself.' + Fore.BLACK)
    
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
        print(Fore.RED + 'Argument D must be of type grumpy dict.' + Fore.BLACK); return
    Type = D.get('Measurement_type', '')
    if not Type == 'fermi_map':
        Type = D.get('Type', '')
        if not Type == 'fermi_map':
            print(Fore.RED + 'The data in the dict is not a Fermi map.' + Fore.BLACK); return
    
    recognized_kwargs = ['figsize', 'cmap', 'linewidth']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    figsize = kwargs.get('figsize', (10,5))
    cmap = kwargs.get('cmap', 'bone_r')
    linewidth = kwargs.get('linewidth', 0.75)

    ENERGY = D.get('x')
    ANGLEX = D.get('z')
    ANGLEY = D.get('y')
    print("\n debug -", ENERGY[0], ANGLEX[0], ANGLEY[0], "\n")     # DEBUG
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
    if not type(D) is dict: print(Fore.RED + 'Argument D must be a grumpy dict.' + Fore.BLACK); return
    #
    recognized_kwargs = ['']
    _kwarg_checker(key_list = recognized_kwargs, **kwargs)

    dim = len(np.shape(D.get('int', np.array([]))))
    if dim < 2:
        print(Fore.RED + 'The data seem to be have too few dimensions. This method requires at least two axes.' + Fore.BLACK); return
    elif dim == 2: 
        _Explore2d(D = D, **kwargs)
    elif dim == 3:
        _Explore3d(D = D, **kwargs)
    else:
        print(Fore.RED + 'The data seem to be have too many dimensions. This method handles maximum three axes.' + Fore.BLACK); return


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

    """
    #
    recognized_kwargs = ['plot', 'shup', 'figsize', 'linewidth', 'linestyle', 'exclude', 'median_filter', 
                         'size', 'mode', 'remove_spikes', 'filter_outliers', 'threshold']
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

    global SHERMAN
    sherman = kwargs.get('sherman', SHERMAN)

    shup = kwargs.get('shup', False)
    #
    if not type(D) is dict:
        print(Fore.RED + 'Argument D must a grumpy dict.' + Fore.BLACK); return {}
    if not D.get('Type', '') == 'spin_edc':
        print(Fore.RED + 'The dicts does not contain spin edc data.' + Fore.BLACK); return D

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
        ax[2].plot(energy, asymmetry, linewidth = linewidth, color = 'k')
        ax[2].set_title('asymmetry')
    
    P1 = (edc_off + edc_on) * (1 + asymmetry / sherman) * 0.5 # there must be a 0.5 here, right?
    P2 = (edc_off + edc_on) * (1 - asymmetry / sherman) * 0.5
    if plot:
        ax[3].plot(energy, P1, color = 'tab:blue', label = 'P1', linewidth = linewidth, linestyle = linestyle)
        ax[3].plot(energy, P2, color = 'tab:red', label = 'P2', linewidth = linewidth, linestyle = linestyle)
        ax[3].legend()
        ax[3].set_title('polarization, S = {0}'.format(sherman))
    
    DD = copy.deepcopy(D)
    DD.update({'int': edc_off})
    DD.update({'int_on': edc_on})
    DD.update({'asymmetry': asymmetry})
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
        print(Fore.RED + 'Argument D must a grumpy dict.' + Fore.BLACK); return D
    if not D.get('Type', '').startswith('spin_deflector'):
        print(Fore.RED + 'The dicts does not contain deflector spin data.' + Fore.BLACK); return D
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
        print(Fore.RED + 'Probably most likely there is potentially some sloppy coding involved here. Sorry. Aborting.' + Fore.BLACK); return {}        

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
        print(Fore.GREEN + 'Please use help(grumpy.Polarization).' + Fore.BLACK)
        return {}
    for i,t in enumerate(types):
        if not t == 'none' and not t in allowed_types: types[i] = 'wrong'
    if any(x=='wrong' for x in types):
        print(Fore.RED + 'One or more of the passed data dicts are of the wrong measurement types for this method.' + Fore.BLACK)
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
        print(Fore.GREEN + 'For just one data set with coil 2 use QuickSpin() or QuickSpinMDC().' + Fore.BLACK)
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
        print(Fore.RED + 'Argument D must be a dict containing fit results.' + Fore.BLACK); return
    if not('yfit' in D and 'par' in D and 'method' in D):
        print(Fore.RED + 'Argument D must be a dict containing fit results.' + Fore.BLACK); return
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
        print(Fore.RED + 'Argument fit must be a dict from e.g. Fit() or similar method.' + Fore.BLACK); return fit
    if not "method" in fit:
        print(Fore.RED + 'Argument fit must contain a reference to the fit profile used.' + Fore.BLACK); return fit
    #
    xaxis = fit.get('x')
    y = fit.get('y')
    if type(xaxis) is type(None):
        print(Fore.RED + 'Argument fit has no key x.' + Fore.BLACK); return fit
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
        print(Fore.RED + 'Argument D must be of type grumpy dict.' + Fore.BLACK); return D
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
        print(Fore.LIGHTRED_EX + "This method is currently only working for 1d data types." + Fore.BLACK)


# ==============================================================================================
# ============================================================================================== S M O O T H
# ==============================================================================================



# ----------------------------------- Remove spikes

# This guy is (e.g.) used when loading spiky spind edc:s.
# Args to pass to Load are filter
def removeSpikes(Y,threshold=2):
    """
    Curtesy of Craig.
    Removes spikes from a 1D intensity array.
    """
    if len(np.shape(Y)) > 2:
        print(Fore.RED + "This works for 1d and 2d data."+ Fore.BLACK)
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


def filterOutliers(Y,threshold=2):
    print(Fore.LIGHTRED_EX + "Note: filterOutliers() is a legacy method name. Will be removed. Use removeSpikes()."+ Fore.BLACK)
    return removeSpikes(Y,threshold=2)



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
            print(Fore.RED + "Argument 'x0' must be a number (or None to be ignored)." + Fore.BLACK); return None, None
        if x0 < x.min() or x0 > x.max():
            print(Fore.RED + f"Argument 'x0' must be in the range ({x.min()}, {x.max()})." + Fore.BLACK); return None, None
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
        print(Fore.RED + 'Argument D must be a grumpy-type dict.' + Fore.BLACK); return D
    if D == {}:
        print(Fore.RED + 'Argument D must be a grumpy-type dict.' + Fore.BLACK); return D
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
        msg = Fore.BLACK + ", ".join(m for m in methods)
        print(Fore.RED + "Argument method must be one of: {0}.".format(msg) + Fore.BLACK)
        return D
    #
    INTKEY = ''
    spin_keys = ['int', 'int_on', 'asymmetry']
    if 'asymmetry' in D: 
        if spin_key == '':
            print("Note: argument spin was not passed so using default spin = 'asymmetry'.")
            spin_key = 'asymmetry'
        if not spin_key in spin_keys:
            msg = Fore.BLACK + ", ".join(m for m in methods)
            print(Fore.RED + "Argument method must be one of: {0}.".format(msg) + Fore.BLACK)
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
        print(Fore.RED + "This method is only accepting 1d and 2d data." + Fore.BLACK)
        return D
    #
    DD = copy.deepcopy(D)
    DD.update({INTKEY: Y})
    return DD

        
    


            
            
            




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
            print(Fore.GREEN + f"{key}" + "{" + Fore.BLACK)
            Dict(item)
            print(Fore.GREEN + '}' + Fore.BLACK)
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
        print(Fore.BLUE + 'Recognized **arguments: {0}'.format(", ".join(k for k in key_list)) + Fore.BLACK)
    c = 0
    for kw in kwargs.keys():
        if not kw in key_list and not kw == 'kw':
            print(Fore.RED + "Argument {0} is not recognized.".format(kw) + Fore.BLACK)
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
            print("for some reason, want it to be like this." + Fore.BLACK)
    if z1[-1] >= z2[0]:
        if not shup: 
            print(Fore.MAGENTA + "Note: is seems like there is an overlap of shift-x axis values. Use the")
            print("SubArray() method to remedy this (if you want)."+ Fore.BLACK)
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

class xps():
    """
    Class containing binding energies for the elements.

    Methods to use:

    Element() or E():           Pass argument 'element' to list binding energies.

    KineticEnergy() or Ek():    Pass arguments 'ek' and 'hv' to list possible elements.
                                Argument 'delta' defines the search window. Argument
                                'order' defines how many orders of hv to search within.
                                Argument filter is a list of elements that limits the search.
    
    BindingEnergy() or Eb():    Pass argument 'eb' to list binding energies within the energy
                                window defined by argument 'delta'. Argument 'filter' is a
                                list of elements that limits the search.
    
    """
    def __init__(self, shup = False):
        if not shup:
            print(Fore.BLUE)
            print("Use .Element() to search for binding energies for an element, .KineticEnergy() to search")
            print("for element binding energies close to a kinetic energy, and .BindingEnergy() to search for")
            print("element binding energies close to a binding energy.")
            print(Fore.MAGENTA + "Credit: All binding energy values comes from Pesto by Craig Polley.")
            print(Fore.BLACK)
        self._make_table()
        self.elements = np.unique(np.transpose(self.xps_table)[1])


    def Element(self, element = 'none'):
        """
        Argument: element as a string, e.g. element = 'Au' lists all binding energies for element Au.
        """
        result = []
        if element in self.elements:
            for i, e in enumerate(self.xps_table):
                if e[1] == element: result.append(self.xps_table[i])
        else:
            if element == 'none':
                print(Fore.RED + "Pass argument element.".format(element) + Fore.BLACK)
            else:
                print(Fore.RED + "Element '{0}' is not in the list.".format(element) + Fore.BLACK)
            return
        for row in result:
            if row[3] == 's': row4 = '   '
            else: row4 = row[4]
            print("{0:<2}   {1:<2}   {2:<1}{3:<1}   {4:<3}  {5:7.1f}".format(row[0], row[1], row[2], row[3], row4, row[5]))

    def KineticEnergy(self, ek = None, hv = None, delta = 5, wf = 4.5, order = 4, filter = []):
        """
        Arguments:  ek as number (kinetic energy), non-optional
                    hv as number (photon energy), non-optional
                    delta as number (energy window) 
                    wf as number (analyzer work function)
                    order as integer (highest order to include)
                    filter as list of strings (filter out elements from the search)
        
        Prints binding energy levels corresponding to the kinetic energy (within a range set by argument delta).
        """
        if type(ek) is type(None) or type(hv) is type(None):
            print(Fore.RED + "Arguments ek and hv must be numbers." + Fore.BLACK)
            return
        result = []
        orders = list(range(1,order+1))
        result_orders = []
        for i, element in enumerate(self.xps_table):
            Eb_tab = float(element[5])
            for order in orders:
                Eb = order * hv - ek - wf
                if Eb_tab - delta <= Eb and Eb <= Eb_tab + delta:
                    append = True
                    if len(filter) > 0:
                        if not element[1] in filter:
                            append = False
                    if append:
                        result.append(self.xps_table[i])
                        result_orders.append(order)
        if len(result) == 0:
            print(Fore.RED + "Could not find any matching results for {0} < Ek < {1}.".format(ek - delta/2, ek + delta/2) + Fore.BLACK)
            return
        for i, row in enumerate(result):
            if row[3] == 's': row4 = '   '
            else: row4 = row[4]
            if result_orders[i] == 1: ord = ''
            else: ord = result_orders[i]
            print("{0:<2}   {1:<2}   {2:<1}{3:<1}   {4:<3}  {5:6.1f}   {6}".format(row[0], row[1], row[2], row[3], row4, row[5], ord))

    def BindingEnergy(self, eb = None, delta = 5, filter = []):
        """
        Arguments:  eb as number (binding energy)
                    delta as number (energy window)
                    filter as list of strings (filter out elements from the search)
        """
        if type(eb) is type(None) or type(delta) is type(None):
            print(Fore.RED + "Arguments eb and delta must be numbers." + Fore.BLACK)
            return
        result = []
        for i, element in enumerate(self.xps_table):
            if (eb - delta/2) <= element[5] and element[5] <= (eb + delta/2):
                append = True
                if len(filter) > 0:
                    if not element[1] in filter:
                        append = False
                if append:
                    result.append(element)
        if len(result) == 0:
            print(Fore.RED + "Could not find any matching results for {0} < Eb < {1}.".format(eb - delta/2, eb + delta/2) + Fore.BLACK)
            return
        for i, row in enumerate(result):
            if row[3] == 's': row4 = '   '
            else: row4 = row[4]
            print("{0:<2}   {1:<2}   {2:<1}{3:<1}   {4:<3}  {5:6.1f}".format(row[0], row[1], row[2], row[3], row4, row[5]))

    def E(self, element = 'none'):
        self.Element(element = element)
    
    def Ek(self, ek = None, hv = None, delta = 5, wf = 4.5, order = 4, filter = []):
        self.KineticEnergy(ek = ek, hv = hv, delta = delta, wf = wf, order = order, filter = filter)
    
    def Eb(self, eb = None, delta = 5, filter = []):
        self.BindingEnergy(self, eb = eb, delta = delta, filter = filter)
    

    def _make_table(self):
        xps_table=[]

        xps_table.append([3, 'Li', '1', 's', '0', 54.7])
        xps_table.append([4, 'Be', '1', 's', '0', 111.5])
        xps_table.append([5, 'B', '1', 's', '0', 188])
        xps_table.append([6, 'C', '1', 's', '0', 284.2])
        xps_table.append([7, 'N', '1', 's', '0', 409.9])
        xps_table.append([7, 'N', '2', 's', '0', 37.3])
        xps_table.append([8, 'O', '1', 's', '0', 543.1])
        xps_table.append([8, 'O', '2', 's', '0', 41.6])
        xps_table.append([9, 'F', '1', 's', '0', 696.7])
        xps_table.append([10, 'Ne', '1', 's', '0', 870.2])
        xps_table.append([10, 'Ne', '2', 's', '0', 48.5])
        xps_table.append([10, 'Ne', '2', 'p', '1/2', 21.7])
        xps_table.append([10, 'Ne', '2', 'p', '3/2', 21.6])


        xps_table.append([11, 'Na', '1', 's', '0', 1070.8])
        xps_table.append([11, 'Na', '2', 's', '0', 63.5])
        xps_table.append([11, 'Na', '2', 'p', '1/2', 30.4])
        xps_table.append([11, 'Na', '2', 'p', '3/2', 30.5])
        xps_table.append([12, 'Mg', '1', 's', '0', 1303.0])
        xps_table.append([12, 'Mg', '2', 's', '0', 88.6])
        xps_table.append([12, 'Mg', '2', 'p', '1/2', 49.6])
        xps_table.append([12, 'Mg', '2', 'p', '3/2', 49.2])
        xps_table.append([13, 'Al', '1', 's', '0', 1559.0])
        xps_table.append([13, 'Al', '2', 's', '0', 117.8])
        xps_table.append([13, 'Al', '2', 'p', '1/2', 72.9])
        xps_table.append([13, 'Al', '2', 'p', '3/2', 72.5])
        xps_table.append([14, 'Si', '1', 's', '0', 1839.0])
        xps_table.append([14, 'Si', '2', 's', '0', 149.7])
        xps_table.append([14, 'Si', '2', 'p', '1/2', 99.8])
        xps_table.append([14, 'Si', '2', 'p', '3/2', 99.2])
        xps_table.append([15, 'P', '1', 's', '0', 2145.5])
        xps_table.append([15, 'P', '2', 's', '0', 189.0])
        xps_table.append([15, 'P', '2', 'p', '1/2', 136.0])
        xps_table.append([15, 'P', '2', 'p', '3/2', 135.0])
        xps_table.append([16, 'S', '1', 's', '0', 2472.0])
        xps_table.append([16, 'S', '2', 's', '0', 230.9])
        xps_table.append([16, 'S', '2', 'p', '1/2', 163.6])
        xps_table.append([16, 'S', '2', 'p', '3/2', 162.5])
        xps_table.append([17, 'Cl', '1', 's', '0', 2822])
        xps_table.append([17, 'Cl', '2', 's', '0', 270])
        xps_table.append([17, 'Cl', '2', 'p', '1/2', 202])
        xps_table.append([17, 'Cl', '2', 'p', '3/2', 200])
        xps_table.append([18, 'Ar', '1', 's', '0', 3205.9])
        xps_table.append([18, 'Ar', '2', 's', '0', 326.3])
        xps_table.append([18, 'Ar', '2', 'p', '1/2', 250.6])
        xps_table.append([18, 'Ar', '2', 'p', '3/2', 248.4])
        xps_table.append([18, 'Ar', '3', 's', '0', 29.3])
        xps_table.append([18, 'Ar', '3', 'p', '1/2', 15.9])
        xps_table.append([18, 'Ar', '3', 'p', '3/2', 15.7])


        xps_table.append([19, 'K', '1', 's', '0', 3608.4])
        xps_table.append([19, 'K', '2', 's', '0', 378.6])
        xps_table.append([19, 'K', '2', 'p', '1/2', 297.3])
        xps_table.append([19, 'K', '2', 'p', '3/2', 294.6])
        xps_table.append([19, 'K', '3', 's', '0', 34.8])
        xps_table.append([19, 'K', '3', 'p', '1/2', 18.3])
        xps_table.append([19, 'K', '3', 'p', '3/2', 18.3])
        xps_table.append([20, 'Ca', '1', 's', '0', 4038.5])
        xps_table.append([20, 'Ca', '2', 's', '0', 438.4])
        xps_table.append([20, 'Ca', '2', 'p', '1/2', 349.7])
        xps_table.append([20, 'Ca', '2', 'p', '3/2', 346.2])
        xps_table.append([20, 'Ca', '3', 's', '0', 44.3])
        xps_table.append([20, 'Ca', '3', 'p', '1/2', 25.4])
        xps_table.append([20, 'Ca', '3', 'p', '3/2', 25.4])
        xps_table.append([21, 'Sc', '1', 's', '0', 4492])
        xps_table.append([21, 'Sc', '2', 's', '0', 498])
        xps_table.append([21, 'Sc', '2', 'p', '1/2', 403.6])
        xps_table.append([21, 'Sc', '2', 'p', '3/2', 398.7])
        xps_table.append([21, 'Sc', '3', 's', '0', 51.1])
        xps_table.append([21, 'Sc', '3', 'p', '1/2', 28.3])
        xps_table.append([21, 'Sc', '3', 'p', '3/2', 28.3])
        xps_table.append([22, 'Ti', '1', 's', '0', 4966])
        xps_table.append([22, 'Ti', '2', 's', '0', 560.9])
        xps_table.append([22, 'Ti', '2', 'p', '1/2', 460.2])
        xps_table.append([22, 'Ti', '2', 'p', '3/2', 453.8])
        xps_table.append([22, 'Ti', '3', 's', '0', 58.7])
        xps_table.append([22, 'Ti', '3', 'p', '1/2', 32.6])
        xps_table.append([22, 'Ti', '3', 'p', '3/2', 32.6])
        xps_table.append([23, 'V', '1', 's', '0', 5465])
        xps_table.append([23, 'V', '2', 's', '0', 626.7])
        xps_table.append([23, 'V', '2', 'p', '1/2', 519.8])
        xps_table.append([23, 'V', '2', 'p', '3/2', 512.1])
        xps_table.append([23, 'V', '3', 's', '0', 66.3])
        xps_table.append([23, 'V', '3', 'p', '1/2', 37.2])
        xps_table.append([23, 'V', '3', 'p', '3/2', 37.2])
        xps_table.append([24, 'Cr', '1', 's', '0', 5989])
        xps_table.append([24, 'Cr', '2', 's', '0', 696])
        xps_table.append([24, 'Cr', '2', 'p', '1/2', 583.8])
        xps_table.append([24, 'Cr', '2', 'p', '3/2', 574.1])
        xps_table.append([24, 'Cr', '3', 's', '0', 74.1])
        xps_table.append([24, 'Cr', '3', 'p', '1/2', 42.2])
        xps_table.append([24, 'Cr', '3', 'p', '3/2', 42.2])
        xps_table.append([25, 'Mn', '1', 's', '0', 6539])
        xps_table.append([25, 'Mn', '2', 's', '0', 769.1])
        xps_table.append([25, 'Mn', '2', 'p', '1/2', 649.9])
        xps_table.append([25, 'Mn', '2', 'p', '3/2', 638.7])
        xps_table.append([25, 'Mn', '3', 's', '0', 82.3])
        xps_table.append([25, 'Mn', '3', 'p', '1/2', 47.2])
        xps_table.append([25, 'Mn', '3', 'p', '3/2', 47.2])
        xps_table.append([26, 'Fe', '1', 's', '0', 7112])
        xps_table.append([26, 'Fe', '2', 's', '0', 844.6])
        xps_table.append([26, 'Fe', '2', 'p', '1/2', 719.9])
        xps_table.append([26, 'Fe', '2', 'p', '3/2', 706.8])
        xps_table.append([26, 'Fe', '3', 's', '0', 91.3])
        xps_table.append([26, 'Fe', '3', 'p', '1/2', 52.7])
        xps_table.append([26, 'Fe', '3', 'p', '3/2', 52.7])
        xps_table.append([27, 'Co', '1', 's', '0', 7709])
        xps_table.append([27, 'Co', '2', 's', '0', 925.1])
        xps_table.append([27, 'Co', '2', 'p', '1/2', 793.2])
        xps_table.append([27, 'Co', '2', 'p', '3/2', 778.1])
        xps_table.append([27, 'Co', '3', 's', '0', 101])
        xps_table.append([27, 'Co', '3', 'p', '1/2', 58.9])
        xps_table.append([27, 'Co', '3', 'p', '3/2', 59.9])
        xps_table.append([28, 'Ni', '1', 's', '0', 8333])
        xps_table.append([28, 'Ni', '2', 's', '0', 1008.6])
        xps_table.append([28, 'Ni', '2', 'p', '1/2', 870])
        xps_table.append([28, 'Ni', '2', 'p', '3/2', 852.7])
        xps_table.append([28, 'Ni', '3', 's', '0', 110.8])
        xps_table.append([28, 'Ni', '3', 'p', '1/2', 68])
        xps_table.append([28, 'Ni', '3', 'p', '3/2', 66.2])
        xps_table.append([29, 'Cu', '1', 's', '0', 8979])
        xps_table.append([29, 'Cu', '2', 's', '0', 1096.7])
        xps_table.append([29, 'Cu', '2', 'p', '1/2', 952.3])
        xps_table.append([29, 'Cu', '2', 'p', '3/2', 932.7])
        xps_table.append([29, 'Cu', '3', 's', '0', 122.5])
        xps_table.append([29, 'Cu', '3', 'p', '1/2', 77.3])
        xps_table.append([29, 'Cu', '3', 'p', '3/2', 75.1])
        xps_table.append([30, 'Zn', '1', 's', '0', 9659])
        xps_table.append([30, 'Zn', '2', 's', '0', 1196.2])
        xps_table.append([30, 'Zn', '2', 'p', '1/2', 1044.9])
        xps_table.append([30, 'Zn', '2', 'p', '3/2', 1021.8])
        xps_table.append([30, 'Zn', '3', 's', '0', 139.8])
        xps_table.append([30, 'Zn', '3', 'p', '1/2', 91.4])
        xps_table.append([30, 'Zn', '3', 'p', '3/2', 88.6])
        xps_table.append([30, 'Zn', '3', 'd', '3/2', 10.2])
        xps_table.append([30, 'Zn', '3', 'd', '5/2', 10.1])
        xps_table.append([31, 'Ga', '1', 's', '0', 10367])
        xps_table.append([31, 'Ga', '2', 's', '0', 1299])
        xps_table.append([31, 'Ga', '2', 'p', '1/2', 1143.2])
        xps_table.append([31, 'Ga', '2', 'p', '3/2', 1116.4])
        xps_table.append([31, 'Ga', '3', 's', '0', 159.5])
        xps_table.append([31, 'Ga', '3', 'p', '1/2', 103.5])
        xps_table.append([31, 'Ga', '3', 'p', '3/2', 100])
        xps_table.append([31, 'Ga', '3', 'd', '3/2', 18.7])
        xps_table.append([31, 'Ga', '3', 'd', '5/2', 18.7])
        xps_table.append([32, 'Ge', '1', 's', '0', 11103])
        xps_table.append([32, 'Ge', '2', 's', '0', 1414.6])
        xps_table.append([32, 'Ge', '2', 'p', '1/2', 1248.1])
        xps_table.append([32, 'Ge', '2', 'p', '3/2', 1217])
        xps_table.append([32, 'Ge', '3', 's', '0', 180.1])
        xps_table.append([32, 'Ge', '3', 'p', '1/2', 124.9])
        xps_table.append([32, 'Ge', '3', 'p', '3/2', 120.8])
        xps_table.append([32, 'Ge', '3', 'd', '3/2', 29.8])
        xps_table.append([32, 'Ge', '3', 'd', '5/2', 29.2])
        xps_table.append([33, 'As', '1', 's', '0', 11867])
        xps_table.append([33, 'As', '2', 's', '0', 1527])
        xps_table.append([33, 'As', '2', 'p', '1/2', 1359.1])
        xps_table.append([33, 'As', '2', 'p', '3/2', 1323.6])
        xps_table.append([33, 'As', '3', 's', '0', 204.7])
        xps_table.append([33, 'As', '3', 'p', '1/2', 146.2])
        xps_table.append([33, 'As', '3', 'p', '3/2', 141.2])
        xps_table.append([33, 'As', '3', 'd', '3/2', 41.7])
        xps_table.append([33, 'As', '3', 'd', '5/2', 41.7])
        xps_table.append([34, 'Se', '1', 's', '0', 12658])
        xps_table.append([34, 'Se', '2', 's', '0', 1652])
        xps_table.append([34, 'Se', '2', 'p', '1/2', 1474.3])
        xps_table.append([34, 'Se', '2', 'p', '3/2', 1433.9])
        xps_table.append([34, 'Se', '3', 's', '0', 229.6])
        xps_table.append([34, 'Se', '3', 'p', '1/2', 166.5])
        xps_table.append([34, 'Se', '3', 'p', '3/2', 160.7])
        xps_table.append([34, 'Se', '3', 'd', '3/2', 55.5])
        xps_table.append([34, 'Se', '3', 'd', '5/2', 54.6])
        xps_table.append([35, 'Br', '1', 's', '0', 13474])
        xps_table.append([35, 'Br', '2', 's', '0', 1782])
        xps_table.append([35, 'Br', '2', 'p', '1/2', 1596])
        xps_table.append([35, 'Br', '2', 'p', '3/2', 1550])
        xps_table.append([35, 'Br', '3', 's', '0', 257])
        xps_table.append([35, 'Br', '3', 'p', '1/2', 189])
        xps_table.append([35, 'Br', '3', 'p', '3/2', 182])
        xps_table.append([35, 'Br', '3', 'd', '3/2', 70])
        xps_table.append([35, 'Br', '3', 'd', '5/2', 69])
        xps_table.append([36, 'Kr', '1', 's', '0', 14326])
        xps_table.append([36, 'Kr', '2', 's', '0', 1921])
        xps_table.append([36, 'Kr', '2', 'p', '1/2', 1730.9])
        xps_table.append([36, 'Kr', '2', 'p', '3/2', 1678.4])
        xps_table.append([36, 'Kr', '3', 's', '0', 292.8])
        xps_table.append([36, 'Kr', '3', 'p', '1/2', 222.2])
        xps_table.append([36, 'Kr', '3', 'p', '3/2', 214.4])
        xps_table.append([36, 'Kr', '3', 'd', '3/2', 95])
        xps_table.append([36, 'Kr', '3', 'd', '5/2', 93.8])
        xps_table.append([36, 'Kr', '4', 's', '0', 27.5])
        xps_table.append([36, 'Kr', '4', 'p', '1/2', 14.1])
        xps_table.append([36, 'Kr', '4', 'p', '3/2', 14.1])


        xps_table.append([37, 'Rb', '1', 's', '0', 15200])
        xps_table.append([37, 'Rb', '2', 's', '0', 2065])
        xps_table.append([37, 'Rb', '2', 'p', '1/2', 1864])
        xps_table.append([37, 'Rb', '2', 'p', '3/2', 1804])
        xps_table.append([37, 'Rb', '3', 's', '0', 326.7])
        xps_table.append([37, 'Rb', '3', 'p', '1/2', 248.7])
        xps_table.append([37, 'Rb', '3', 'p', '3/2', 239.1])
        xps_table.append([37, 'Rb', '3', 'd', '3/2', 113])
        xps_table.append([37, 'Rb', '3', 'd', '5/2', 112])
        xps_table.append([37, 'Rb', '4', 's', '0', 30.5])
        xps_table.append([37, 'Rb', '4', 'p', '1/2', 16.3])
        xps_table.append([37, 'Rb', '4', 'p', '3/2', 15.3])

        xps_table.append([38, 'Sr', '1', 's', '0', 16105])
        xps_table.append([38, 'Sr', '2', 's', '0', 2216])
        xps_table.append([38, 'Sr', '2', 'p', '1/2', 2007])
        xps_table.append([38, 'Sr', '2', 'p', '3/2', 1940])
        xps_table.append([38, 'Sr', '3', 's', '0', 358.7])
        xps_table.append([38, 'Sr', '3', 'p', '1/2', 280.3])
        xps_table.append([38, 'Sr', '3', 'p', '3/2', 270])
        xps_table.append([38, 'Sr', '3', 'd', '3/2', 136])
        xps_table.append([38, 'Sr', '3', 'd', '5/2', 134.2])
        xps_table.append([38, 'Sr', '4', 's', '0', 38.9])
        xps_table.append([38, 'Sr', '4', 'p', '1/2', 21.6])
        xps_table.append([38, 'Sr', '4', 'p', '3/2', 20.1])

        xps_table.append([39, 'Y', '1', 's', '0', 17038])
        xps_table.append([39, 'Y', '2', 's', '0', 2373])
        xps_table.append([39, 'Y', '2', 'p', '1/2', 2156])
        xps_table.append([39, 'Y', '2', 'p', '3/2', 2080])
        xps_table.append([39, 'Y', '3', 's', '0', 392])
        xps_table.append([39, 'Y', '3', 'p', '1/2', 310.6])
        xps_table.append([39, 'Y', '3', 'p', '3/2', 298.8])
        xps_table.append([39, 'Y', '3', 'd', '3/2', 157.7])
        xps_table.append([39, 'Y', '3', 'd', '5/2', 155.8])
        xps_table.append([39, 'Y', '4', 's', '0', 43.8])
        xps_table.append([39, 'Y', '4', 'p', '1/2', 24.4])
        xps_table.append([39, 'Y', '4', 'p', '3/2', 23.1])

        xps_table.append([40, 'Zr', '1', 's', '0', 17998])
        xps_table.append([40, 'Zr', '2', 's', '0', 2532])
        xps_table.append([40, 'Zr', '2', 'p', '1/2', 2307])
        xps_table.append([40, 'Zr', '2', 'p', '3/2', 2223])
        xps_table.append([40, 'Zr', '3', 's', '0', 430.3])
        xps_table.append([40, 'Zr', '3', 'p', '1/2', 343.5])
        xps_table.append([40, 'Zr', '3', 'p', '3/2', 329.8])
        xps_table.append([40, 'Zr', '3', 'd', '3/2', 181.1])
        xps_table.append([40, 'Zr', '3', 'd', '5/2', 178.8])
        xps_table.append([40, 'Zr', '4', 's', '0', 50.6])
        xps_table.append([40, 'Zr', '4', 'p', '1/2', 28.5])
        xps_table.append([40, 'Zr', '4', 'p', '3/2', 27.1])

        xps_table.append([41, 'Nb', '1', 's', '0', 18986])
        xps_table.append([41, 'Nb', '2', 's', '0', 2698])
        xps_table.append([41, 'Nb', '2', 'p', '1/2', 2465])
        xps_table.append([41, 'Nb', '2', 'p', '3/2', 2371])
        xps_table.append([41, 'Nb', '3', 's', '0', 466.6])
        xps_table.append([41, 'Nb', '3', 'p', '1/2', 376.1])
        xps_table.append([41, 'Nb', '3', 'p', '3/2', 360.6])
        xps_table.append([41, 'Nb', '3', 'd', '3/2', 205])
        xps_table.append([41, 'Nb', '3', 'd', '5/2', 202.3])
        xps_table.append([41, 'Nb', '4', 's', '0', 56.4])
        xps_table.append([41, 'Nb', '4', 'p', '1/2', 32.6])
        xps_table.append([41, 'Nb', '4', 'p', '3/2', 30.8])

        xps_table.append([42, 'Mo', '1', 's', '0', 20000])
        xps_table.append([42, 'Mo', '2', 's', '0', 2866])
        xps_table.append([42, 'Mo', '2', 'p', '1/2', 2625])
        xps_table.append([42, 'Mo', '2', 'p', '3/2', 2520])
        xps_table.append([42, 'Mo', '3', 's', '0', 506.3])
        xps_table.append([42, 'Mo', '3', 'p', '1/2', 411.6])
        xps_table.append([42, 'Mo', '3', 'p', '3/2', 394])
        xps_table.append([42, 'Mo', '3', 'd', '3/2', 231.1])
        xps_table.append([42, 'Mo', '3', 'd', '5/2', 227.9])
        xps_table.append([42, 'Mo', '4', 's', '0', 63.2])
        xps_table.append([42, 'Mo', '4', 'p', '1/2', 37.6])
        xps_table.append([42, 'Mo', '4', 'p', '3/2', 35.5])

        xps_table.append([43, 'Tc', '1', 's', '0', 21044])
        xps_table.append([43, 'Tc', '2', 's', '0', 3043])
        xps_table.append([43, 'Tc', '2', 'p', '1/2', 2793])
        xps_table.append([43, 'Tc', '2', 'p', '3/2', 2677])
        xps_table.append([43, 'Tc', '3', 's', '0', 544])
        xps_table.append([43, 'Tc', '3', 'p', '1/2', 447.6])
        xps_table.append([43, 'Tc', '3', 'p', '3/2', 417.7])
        xps_table.append([43, 'Tc', '3', 'd', '3/2', 257.6])
        xps_table.append([43, 'Tc', '3', 'd', '5/2', 253.9])
        xps_table.append([43, 'Tc', '4', 's', '0', 69.5])
        xps_table.append([43, 'Tc', '4', 'p', '1/2', 42.3])
        xps_table.append([43, 'Tc', '4', 'p', '3/2', 39.9])

        xps_table.append([44, 'Ru', '1', 's', '0', 22117])
        xps_table.append([44, 'Ru', '2', 's', '0', 3224])
        xps_table.append([44, 'Ru', '2', 'p', '1/2', 2967])
        xps_table.append([44, 'Ru', '2', 'p', '3/2', 2838])
        xps_table.append([44, 'Ru', '3', 's', '0', 586.1])
        xps_table.append([44, 'Ru', '3', 'p', '1/2', 483.3])
        xps_table.append([44, 'Ru', '3', 'p', '3/2', 461.5])
        xps_table.append([44, 'Ru', '3', 'd', '3/2', 284.2])
        xps_table.append([44, 'Ru', '3', 'd', '5/2', 280])
        xps_table.append([44, 'Ru', '4', 's', '0', 75])
        xps_table.append([44, 'Ru', '4', 'p', '1/2', 46.3])
        xps_table.append([44, 'Ru', '4', 'p', '3/2', 43.2])

        xps_table.append([45, 'Rh', '1', 's', '0', 23220])
        xps_table.append([45, 'Rh', '2', 's', '0', 3412])
        xps_table.append([45, 'Rh', '2', 'p', '1/2', 3146])
        xps_table.append([45, 'Rh', '2', 'p', '3/2', 3004])
        xps_table.append([45, 'Rh', '3', 's', '0', 628.1])
        xps_table.append([45, 'Rh', '3', 'p', '1/2', 521.3])
        xps_table.append([45, 'Rh', '3', 'p', '3/2', 496.5])
        xps_table.append([45, 'Rh', '3', 'd', '3/2', 311.9])
        xps_table.append([45, 'Rh', '3', 'd', '5/2', 307.2])
        xps_table.append([45, 'Rh', '4', 's', '0', 81.4])
        xps_table.append([45, 'Rh', '4', 'p', '1/2', 50.5])
        xps_table.append([45, 'Rh', '4', 'p', '3/2', 47.3])

        xps_table.append([46, 'Pd', '1', 's', '0', 24350])
        xps_table.append([46, 'Pd', '2', 's', '0', 3604])
        xps_table.append([46, 'Pd', '2', 'p', '1/2', 3330])
        xps_table.append([46, 'Pd', '2', 'p', '3/2', 3173])
        xps_table.append([46, 'Pd', '3', 's', '0', 671.6])
        xps_table.append([46, 'Pd', '3', 'p', '1/2', 559.9])
        xps_table.append([46, 'Pd', '3', 'p', '3/2', 532.3])
        xps_table.append([46, 'Pd', '3', 'd', '3/2', 340.5])
        xps_table.append([46, 'Pd', '3', 'd', '5/2', 335.2])
        xps_table.append([46, 'Pd', '4', 's', '0', 87.1])
        xps_table.append([46, 'Pd', '4', 'p', '1/2', 55.7])
        xps_table.append([46, 'Pd', '4', 'p', '3/2', 50.9])

        xps_table.append([47, 'Ag', '1', 's', '0', 25514])
        xps_table.append([47, 'Ag', '2', 's', '0', 3806])
        xps_table.append([47, 'Ag', '2', 'p', '1/2', 3524])
        xps_table.append([47, 'Ag', '2', 'p', '3/2', 3351])
        xps_table.append([47, 'Ag', '3', 's', '0', 719])
        xps_table.append([47, 'Ag', '3', 'p', '1/2', 603.8])
        xps_table.append([47, 'Ag', '3', 'p', '3/2', 573])
        xps_table.append([47, 'Ag', '3', 'd', '3/2', 374])
        xps_table.append([47, 'Ag', '3', 'd', '5/2', 368.3])
        xps_table.append([47, 'Ag', '4', 's', '0', 97])
        xps_table.append([47, 'Ag', '4', 'p', '1/2', 63.7])
        xps_table.append([47, 'Ag', '4', 'p', '3/2', 58.3])

        xps_table.append([48, 'Cd', '1', 's', '0', 26711])
        xps_table.append([48, 'Cd', '2', 's', '0', 4018])
        xps_table.append([48, 'Cd', '2', 'p', '1/2', 3727])
        xps_table.append([48, 'Cd', '2', 'p', '3/2', 3538])
        xps_table.append([48, 'Cd', '3', 's', '0', 772])
        xps_table.append([48, 'Cd', '3', 'p', '1/2', 652.6])
        xps_table.append([48, 'Cd', '3', 'p', '3/2', 618.4])
        xps_table.append([48, 'Cd', '3', 'd', '3/2', 411.9])
        xps_table.append([48, 'Cd', '3', 'd', '5/2', 405.2])
        xps_table.append([48, 'Cd', '4', 's', '0', 109.8])
        xps_table.append([48, 'Cd', '4', 'p', '1/2', 63.9])
        xps_table.append([48, 'Cd', '4', 'p', '3/2', 63.9])
        xps_table.append([48, 'Cd', '4', 'd', '3/2', 11.7])
        xps_table.append([48, 'Cd', '4', 'd', '5/2', 10.7])

        xps_table.append([49, 'In', '1', 's', '0', 27940])
        xps_table.append([49, 'In', '2', 's', '0', 4238])
        xps_table.append([49, 'In', '2', 'p', '1/2', 3938])
        xps_table.append([49, 'In', '2', 'p', '3/2', 3730])
        xps_table.append([49, 'In', '3', 's', '0', 827.2])
        xps_table.append([49, 'In', '3', 'p', '1/2', 703.2])
        xps_table.append([49, 'In', '3', 'p', '3/2', 665.3])
        xps_table.append([49, 'In', '3', 'd', '3/2', 451.4])
        xps_table.append([49, 'In', '3', 'd', '5/2', 443.9])
        xps_table.append([49, 'In', '4', 's', '0', 122.9])
        xps_table.append([49, 'In', '4', 'p', '1/2', 73.5])
        xps_table.append([49, 'In', '4', 'p', '3/2', 73.5])
        xps_table.append([49, 'In', '4', 'd', '3/2', 17.7])
        xps_table.append([49, 'In', '4', 'd', '5/2', 16.9])

        xps_table.append([50, 'Sn', '1', 's', '0', 29200])
        xps_table.append([50, 'Sn', '2', 's', '0', 4465])
        xps_table.append([50, 'Sn', '2', 'p', '1/2', 4156])
        xps_table.append([50, 'Sn', '2', 'p', '3/2', 3929])
        xps_table.append([50, 'Sn', '3', 's', '0', 884.7])
        xps_table.append([50, 'Sn', '3', 'p', '1/2', 756.5])
        xps_table.append([50, 'Sn', '3', 'p', '3/2', 714.6])
        xps_table.append([50, 'Sn', '3', 'd', '3/2', 493.2])
        xps_table.append([50, 'Sn', '3', 'd', '5/2', 484.9])
        xps_table.append([50, 'Sn', '4', 's', '0', 137.1])
        xps_table.append([50, 'Sn', '4', 'p', '1/2', 83.6])
        xps_table.append([50, 'Sn', '4', 'p', '3/2', 83.6])
        xps_table.append([50, 'Sn', '4', 'd', '3/2', 24.9])
        xps_table.append([50, 'Sn', '4', 'd', '5/2', 23.9])

        xps_table.append([51, 'Sb', '1', 's', '0', 30491])
        xps_table.append([51, 'Sb', '2', 's', '0', 4698])
        xps_table.append([51, 'Sb', '2', 'p', '1/2', 4380])
        xps_table.append([51, 'Sb', '2', 'p', '3/2', 4132])
        xps_table.append([51, 'Sb', '3', 's', '0', 940])
        xps_table.append([51, 'Sb', '3', 'p', '1/2', 812.7])
        xps_table.append([51, 'Sb', '3', 'p', '3/2', 766.4])
        xps_table.append([51, 'Sb', '3', 'd', '3/2', 537.5])
        xps_table.append([51, 'Sb', '3', 'd', '5/2', 528.2])
        xps_table.append([51, 'Sb', '4', 's', '0', 153.2])
        xps_table.append([51, 'Sb', '4', 'p', '1/2', 95.6])
        xps_table.append([51, 'Sb', '4', 'p', '3/2', 95.6])
        xps_table.append([51, 'Sb', '4', 'd', '3/2', 33.3])
        xps_table.append([51, 'Sb', '4', 'd', '5/2', 32.1])

        xps_table.append([52, 'Te', '1', 's', '0', 31814])
        xps_table.append([52, 'Te', '2', 's', '0', 4939])
        xps_table.append([52, 'Te', '2', 'p', '1/2', 4612])
        xps_table.append([52, 'Te', '2', 'p', '3/2', 4341])
        xps_table.append([52, 'Te', '3', 's', '0', 1006])
        xps_table.append([52, 'Te', '3', 'p', '1/2', 870.8])
        xps_table.append([52, 'Te', '3', 'p', '3/2', 820.8])
        xps_table.append([52, 'Te', '3', 'd', '3/2', 583.4])
        xps_table.append([52, 'Te', '3', 'd', '5/2', 573])
        xps_table.append([52, 'Te', '4', 's', '0', 169.4])
        xps_table.append([52, 'Te', '4', 'p', '1/2', 103.3])
        xps_table.append([52, 'Te', '4', 'p', '3/2', 103.3])
        xps_table.append([52, 'Te', '4', 'd', '3/2', 41.9])
        xps_table.append([52, 'Te', '4', 'd', '5/2', 40.4])

        xps_table.append([53, 'I', '1', 's', '0', 33169])
        xps_table.append([53, 'I', '2', 's', '0', 5188])
        xps_table.append([53, 'I', '2', 'p', '1/2', 4852])
        xps_table.append([53, 'I', '2', 'p', '3/2', 4557])
        xps_table.append([53, 'I', '3', 's', '0', 1072])
        xps_table.append([53, 'I', '3', 'p', '1/2', 931])
        xps_table.append([53, 'I', '3', 'p', '3/2', 875])
        xps_table.append([53, 'I', '3', 'd', '3/2', 630.8])
        xps_table.append([53, 'I', '3', 'd', '5/2', 619.3])
        xps_table.append([53, 'I', '4', 's', '0', 186])
        xps_table.append([53, 'I', '4', 'p', '1/2', 123])
        xps_table.append([53, 'I', '4', 'p', '3/2', 123])
        xps_table.append([53, 'I', '4', 'd', '3/2', 50.6])
        xps_table.append([53, 'I', '4', 'd', '5/2', 48.9])

        xps_table.append([54, 'Xe', '1', 's', '0', 34561])
        xps_table.append([54, 'Xe', '2', 's', '0', 5453])
        xps_table.append([54, 'Xe', '2', 'p', '1/2', 5107])
        xps_table.append([54, 'Xe', '2', 'p', '3/2', 4786])
        xps_table.append([54, 'Xe', '3', 's', '0', 1148.7])
        xps_table.append([54, 'Xe', '3', 'p', '1/2', 1002.1])
        xps_table.append([54, 'Xe', '3', 'p', '3/2', 940.6])
        xps_table.append([54, 'Xe', '3', 'd', '3/2', 689])
        xps_table.append([54, 'Xe', '3', 'd', '5/2', 676.4])
        xps_table.append([54, 'Xe', '4', 's', '0', 213.2])
        xps_table.append([54, 'Xe', '4', 'p', '1/2', 146.7])
        xps_table.append([54, 'Xe', '4', 'p', '3/2', 145.5])
        xps_table.append([54, 'Xe', '4', 'd', '3/2', 69.5])
        xps_table.append([54, 'Xe', '4', 'd', '5/2', 67.5])
        xps_table.append([54, 'Xe', '5', 's', '0', 23.3])
        xps_table.append([54, 'Xe', '5', 'p', '1/2', 13.4])
        xps_table.append([54, 'Xe', '5', 'p', '3/2', 12.1])


        xps_table.append([55, 'Cs', '1', 's', '0', 35985])
        xps_table.append([55, 'Cs', '2', 's', '0', 5714])
        xps_table.append([55, 'Cs', '2', 'p', '1/2', 5359])
        xps_table.append([55, 'Cs', '2', 'p', '3/2', 5012])
        xps_table.append([55, 'Cs', '3', 's', '0', 1211])
        xps_table.append([55, 'Cs', '3', 'p', '1/2', 1071])
        xps_table.append([55, 'Cs', '3', 'p', '3/2', 1003])
        xps_table.append([55, 'Cs', '3', 'd', '3/2', 740.5])
        xps_table.append([55, 'Cs', '3', 'd', '5/2', 726.6])
        xps_table.append([55, 'Cs', '4', 's', '0', 232.3])
        xps_table.append([55, 'Cs', '4', 'p', '1/2', 172.4])
        xps_table.append([55, 'Cs', '4', 'p', '3/2', 161.3])
        xps_table.append([55, 'Cs', '4', 'd', '3/2', 79.8])
        xps_table.append([55, 'Cs', '4', 'd', '5/2', 77.5])
        xps_table.append([55, 'Cs', '5', 's', '0', 22.7])
        xps_table.append([55, 'Cs', '5', 'p', '1/2', 14.2])
        xps_table.append([55, 'Cs', '5', 'p', '3/2', 12.1])

        xps_table.append([56, 'Ba', '1', 's', '0', 37441])
        xps_table.append([56, 'Ba', '2', 's', '0', 5989])
        xps_table.append([56, 'Ba', '2', 'p', '1/2', 5624])
        xps_table.append([56, 'Ba', '2', 'p', '3/2', 5247])
        xps_table.append([56, 'Ba', '3', 's', '0', 1293])
        xps_table.append([56, 'Ba', '3', 'p', '1/2', 1137])
        xps_table.append([56, 'Ba', '3', 'p', '3/2', 1063])
        xps_table.append([56, 'Ba', '3', 'd', '3/2', 795.7])
        xps_table.append([56, 'Ba', '3', 'd', '5/2', 780.5])
        xps_table.append([56, 'Ba', '4', 's', '0', 253.5])
        xps_table.append([56, 'Ba', '4', 'p', '1/2', 192])
        xps_table.append([56, 'Ba', '4', 'p', '3/2', 178.6])
        xps_table.append([56, 'Ba', '4', 'd', '3/2', 92.6])
        xps_table.append([56, 'Ba', '4', 'd', '5/2', 89.9])
        xps_table.append([56, 'Ba', '5', 's', '0', 30.3])
        xps_table.append([56, 'Ba', '5', 'p', '1/2', 17])
        xps_table.append([56, 'Ba', '5', 'p', '3/2', 14.8])

        xps_table.append([57, 'La', '1', 's', '0', 38925])
        xps_table.append([57, 'La', '2', 's', '0', 6266])
        xps_table.append([57, 'La', '2', 'p', '1/2', 5891])
        xps_table.append([57, 'La', '2', 'p', '3/2', 5483])
        xps_table.append([57, 'La', '3', 's', '0', 1362])
        xps_table.append([57, 'La', '3', 'p', '1/2', 1209])
        xps_table.append([57, 'La', '3', 'p', '3/2', 1128])
        xps_table.append([57, 'La', '3', 'd', '3/2', 853])
        xps_table.append([57, 'La', '3', 'd', '5/2', 836])
        xps_table.append([57, 'La', '4', 's', '0', 274.7])
        xps_table.append([57, 'La', '4', 'p', '1/2', 205.8])
        xps_table.append([57, 'La', '4', 'p', '3/2', 196])
        xps_table.append([57, 'La', '4', 'd', '3/2', 105.3])
        xps_table.append([57, 'La', '4', 'd', '5/2', 102.5])
        xps_table.append([57, 'La', '5', 's', '0', 34.3])
        xps_table.append([57, 'La', '5', 'p', '1/2', 19.3])
        xps_table.append([57, 'La', '5', 'p', '3/2', 16.8])

        xps_table.append([58, 'Ce', '1', 's', '0', 40443])
        xps_table.append([58, 'Ce', '2', 's', '0', 6548])
        xps_table.append([58, 'Ce', '2', 'p', '1/2', 6164])
        xps_table.append([58, 'Ce', '2', 'p', '3/2', 5723])
        xps_table.append([58, 'Ce', '3', 's', '0', 1436])
        xps_table.append([58, 'Ce', '3', 'p', '1/2', 1274])
        xps_table.append([58, 'Ce', '3', 'p', '3/2', 1187])
        xps_table.append([58, 'Ce', '3', 'd', '3/2', 902.4])
        xps_table.append([58, 'Ce', '3', 'd', '5/2', 883.8])
        xps_table.append([58, 'Ce', '4', 's', '0', 291])
        xps_table.append([58, 'Ce', '4', 'p', '1/2', 223.2])
        xps_table.append([58, 'Ce', '4', 'p', '3/2', 206.5])
        xps_table.append([58, 'Ce', '4', 'd', '3/2', 109])
        xps_table.append([58, 'Ce', '4', 'f', '5/2', 0.1])
        xps_table.append([58, 'Ce', '4', 'f', '7/2', 0.1])
        xps_table.append([58, 'Ce', '5', 's', '0', 37.8])
        xps_table.append([58, 'Ce', '5', 'p', '1/2', 19.8])
        xps_table.append([58, 'Ce', '5', 'p', '3/2', 17])

        xps_table.append([59, 'Pr', '1', 's', '0', 41991])
        xps_table.append([59, 'Pr', '2', 's', '0', 6835])
        xps_table.append([59, 'Pr', '2', 'p', '1/2', 6440])
        xps_table.append([59, 'Pr', '2', 'p', '3/2', 5964])
        xps_table.append([59, 'Pr', '3', 's', '0', 1511])
        xps_table.append([59, 'Pr', '3', 'p', '1/2', 1337])
        xps_table.append([59, 'Pr', '3', 'p', '3/2', 1242])
        xps_table.append([59, 'Pr', '3', 'd', '3/2', 948.3])
        xps_table.append([59, 'Pr', '3', 'd', '5/2', 928.8])
        xps_table.append([59, 'Pr', '4', 's', '0', 304.5])
        xps_table.append([59, 'Pr', '4', 'p', '1/2', 236.3])
        xps_table.append([59, 'Pr', '4', 'p', '3/2', 217.6])
        xps_table.append([59, 'Pr', '4', 'd', '3/2', 115.1])
        xps_table.append([59, 'Pr', '4', 'd', '5/2', 115.1])
        xps_table.append([59, 'Pr', '4', 'f', '5/2', 2])
        xps_table.append([59, 'Pr', '4', 'f', '7/2', 2])
        xps_table.append([59, 'Pr', '5', 's', '0', 37.4])
        xps_table.append([59, 'Pr', '5', 'p', '1/2', 22.3])
        xps_table.append([59, 'Pr', '5', 'p', '3/2', 22.3])

        xps_table.append([60, 'Nd', '1', 's', '0', 43569])
        xps_table.append([60, 'Nd', '2', 's', '0', 7126])
        xps_table.append([60, 'Nd', '2', 'p', '1/2', 6722])
        xps_table.append([60, 'Nd', '2', 'p', '3/2', 6208])
        xps_table.append([60, 'Nd', '3', 's', '0', 1575])
        xps_table.append([60, 'Nd', '3', 'p', '1/2', 1403])
        xps_table.append([60, 'Nd', '3', 'p', '3/2', 1297])
        xps_table.append([60, 'Nd', '3', 'd', '3/2', 1003.3])
        xps_table.append([60, 'Nd', '3', 'd', '5/2', 980.4])
        xps_table.append([60, 'Nd', '4', 's', '0', 319.2])
        xps_table.append([60, 'Nd', '4', 'p', '1/2', 243.3])
        xps_table.append([60, 'Nd', '4', 'p', '3/2', 224.6])
        xps_table.append([60, 'Nd', '4', 'd', '3/2', 120.5])
        xps_table.append([60, 'Nd', '4', 'd', '5/2', 120.5])
        xps_table.append([60, 'Nd', '4', 'f', '5/2', 1.5])
        xps_table.append([60, 'Nd', '4', 'f', '7/2', 1.5])
        xps_table.append([60, 'Nd', '5', 's', '0', 37.5])
        xps_table.append([60, 'Nd', '5', 'p', '1/2', 21.1])
        xps_table.append([60, 'Nd', '5', 'p', '3/2', 21.1])

        xps_table.append([61, 'Pm', '1', 's', '0', 45184])
        xps_table.append([61, 'Pm', '2', 's', '0', 7428])
        xps_table.append([61, 'Pm', '2', 'p', '1/2', 7013])
        xps_table.append([61, 'Pm', '2', 'p', '3/2', 6459])
        xps_table.append([61, 'Pm', '3', 'p', '1/2', 1471.4])
        xps_table.append([61, 'Pm', '3', 'p', '3/2', 1357])
        xps_table.append([61, 'Pm', '3', 'd', '3/2', 1052])
        xps_table.append([61, 'Pm', '3', 'd', '5/2', 1027])
        xps_table.append([61, 'Pm', '4', 'p', '1/2', 242])
        xps_table.append([61, 'Pm', '4', 'p', '3/2', 242])
        xps_table.append([61, 'Pm', '4', 'd', '3/2', 120])
        xps_table.append([61, 'Pm', '4', 'd', '5/2', 120])

        xps_table.append([62, 'Sm', '1', 's', '0', 46834])
        xps_table.append([62, 'Sm', '2', 's', '0', 7737])
        xps_table.append([62, 'Sm', '2', 'p', '1/2', 7312])
        xps_table.append([62, 'Sm', '2', 'p', '3/2', 6716])
        xps_table.append([62, 'Sm', '3', 's', '0', 1723])
        xps_table.append([62, 'Sm', '3', 'p', '1/2', 1541])
        xps_table.append([62, 'Sm', '3', 'p', '3/2', 1419.8])
        xps_table.append([62, 'Sm', '3', 'd', '3/2', 1110.9])
        xps_table.append([62, 'Sm', '3', 'd', '5/2', 1083.4])
        xps_table.append([62, 'Sm', '4', 's', '0', 347.2])
        xps_table.append([62, 'Sm', '4', 'p', '1/2', 265.6])
        xps_table.append([62, 'Sm', '4', 'p', '3/2', 247.4])
        xps_table.append([62, 'Sm', '4', 'd', '3/2', 129])
        xps_table.append([62, 'Sm', '4', 'd', '5/2', 129])
        xps_table.append([62, 'Sm', '4', 'f', '5/2', 5.2])
        xps_table.append([62, 'Sm', '4', 'f', '7/2', 5.2])
        xps_table.append([62, 'Sm', '5', 's', '0', 37.4])
        xps_table.append([62, 'Sm', '5', 'p', '1/2', 21.3])
        xps_table.append([62, 'Sm', '5', 'p', '3/2', 21.3])

        xps_table.append([63, 'Eu', '1', 's', '0', 48519])
        xps_table.append([63, 'Eu', '2', 's', '0', 8052])
        xps_table.append([63, 'Eu', '2', 'p', '1/2', 7617])
        xps_table.append([63, 'Eu', '2', 'p', '3/2', 6977])
        xps_table.append([63, 'Eu', '3', 's', '0', 1800])
        xps_table.append([63, 'Eu', '3', 'p', '1/2', 1614])
        xps_table.append([63, 'Eu', '3', 'p', '3/2', 1481])
        xps_table.append([63, 'Eu', '3', 'd', '3/2', 1158.6])
        xps_table.append([63, 'Eu', '3', 'd', '5/2', 1127.5])
        xps_table.append([63, 'Eu', '4', 's', '0', 360])
        xps_table.append([63, 'Eu', '4', 'p', '1/2', 284])
        xps_table.append([63, 'Eu', '4', 'p', '3/2', 257])
        xps_table.append([63, 'Eu', '4', 'd', '3/2', 133])
        xps_table.append([63, 'Eu', '4', 'd', '5/2', 127.7])
        xps_table.append([63, 'Eu', '4', 'f', '5/2', 0])
        xps_table.append([63, 'Eu', '4', 'f', '7/2', 0])
        xps_table.append([63, 'Eu', '5', 's', '0', 32])
        xps_table.append([63, 'Eu', '5', 'p', '1/2', 22])
        xps_table.append([63, 'Eu', '5', 'p', '3/2', 22])

        xps_table.append([64, 'Gd', '1', 's', '0', 50239])
        xps_table.append([64, 'Gd', '2', 's', '0', 8376])
        xps_table.append([64, 'Gd', '2', 'p', '1/2', 7930])
        xps_table.append([64, 'Gd', '2', 'p', '3/2', 7243])
        xps_table.append([64, 'Gd', '3', 's', '0', 1881])
        xps_table.append([64, 'Gd', '3', 'p', '1/2', 1688])
        xps_table.append([64, 'Gd', '3', 'p', '3/2', 1544])
        xps_table.append([64, 'Gd', '3', 'd', '3/2', 1221.9])
        xps_table.append([64, 'Gd', '3', 'd', '5/2', 1189.6])
        xps_table.append([64, 'Gd', '4', 's', '0', 378.6])
        xps_table.append([64, 'Gd', '4', 'p', '1/2', 286])
        xps_table.append([64, 'Gd', '4', 'p', '3/2', 271])
        xps_table.append([64, 'Gd', '4', 'd', '5/2', 142.6])
        xps_table.append([64, 'Gd', '4', 'f', '5/2', 8.6])
        xps_table.append([64, 'Gd', '4', 'f', '7/2', 8.6])
        xps_table.append([64, 'Gd', '5', 's', '0', 36])
        xps_table.append([64, 'Gd', '5', 'p', '1/2', 20])
        xps_table.append([64, 'Gd', '5', 'p', '3/2', 20])

        xps_table.append([65, 'Tb', '1', 's', '0', 51996])
        xps_table.append([65, 'Tb', '2', 's', '0', 8708])
        xps_table.append([65, 'Tb', '2', 'p', '1/2', 8252])
        xps_table.append([65, 'Tb', '2', 'p', '3/2', 7514])
        xps_table.append([65, 'Tb', '3', 's', '0', 1968])
        xps_table.append([65, 'Tb', '3', 'p', '1/2', 1768])
        xps_table.append([65, 'Tb', '3', 'p', '3/2', 1611])
        xps_table.append([65, 'Tb', '3', 'd', '3/2', 1276.9])
        xps_table.append([65, 'Tb', '3', 'd', '5/2', 1241.1])
        xps_table.append([65, 'Tb', '4', 's', '0', 396])
        xps_table.append([65, 'Tb', '4', 'p', '1/2', 322.4])
        xps_table.append([65, 'Tb', '4', 'p', '3/2', 284.1])
        xps_table.append([65, 'Tb', '4', 'd', '3/2', 150.5])
        xps_table.append([65, 'Tb', '4', 'd', '5/2', 150.5])
        xps_table.append([65, 'Tb', '4', 'f', '5/2', 7.7])
        xps_table.append([65, 'Tb', '4', 'f', '7/2', 2.4])
        xps_table.append([65, 'Tb', '5', 's', '0', 45.6])
        xps_table.append([65, 'Tb', '5', 'p', '1/2', 28.7])
        xps_table.append([65, 'Tb', '5', 'p', '3/2', 22.6])

        xps_table.append([66, 'Dy', '1', 's', '0', 53789])
        xps_table.append([66, 'Dy', '2', 's', '0', 9046])
        xps_table.append([66, 'Dy', '2', 'p', '1/2', 8581])
        xps_table.append([66, 'Dy', '2', 'p', '3/2', 7790])
        xps_table.append([66, 'Dy', '3', 's', '0', 2047])
        xps_table.append([66, 'Dy', '3', 'p', '1/2', 1842])
        xps_table.append([66, 'Dy', '3', 'p', '3/2', 1676])
        xps_table.append([66, 'Dy', '3', 'd', '3/2', 1333])
        xps_table.append([66, 'Dy', '3', 'd', '5/2', 1292])
        xps_table.append([66, 'Dy', '4', 's', '0', 414.2])
        xps_table.append([66, 'Dy', '4', 'p', '1/2', 333.5])
        xps_table.append([66, 'Dy', '4', 'p', '3/2', 293.2])
        xps_table.append([66, 'Dy', '4', 'd', '3/2', 153.6])
        xps_table.append([66, 'Dy', '4', 'd', '5/2', 153.6])
        xps_table.append([66, 'Dy', '4', 'f', '5/2', 8])
        xps_table.append([66, 'Dy', '4', 'f', '7/2', 4.3])
        xps_table.append([66, 'Dy', '5', 's', '0', 49.9])
        xps_table.append([66, 'Dy', '5', 'p', '1/2', 26.3])
        xps_table.append([66, 'Dy', '5', 'p', '3/2', 26.3])

        xps_table.append([67, 'Ho', '1', 's', '0', 55618])
        xps_table.append([67, 'Ho', '2', 's', '0', 9394])
        xps_table.append([67, 'Ho', '2', 'p', '1/2', 8918])
        xps_table.append([67, 'Ho', '2', 'p', '3/2', 8071])
        xps_table.append([67, 'Ho', '3', 's', '0', 2128])
        xps_table.append([67, 'Ho', '3', 'p', '1/2', 1923])
        xps_table.append([67, 'Ho', '3', 'p', '3/2', 1741])
        xps_table.append([67, 'Ho', '3', 'd', '3/2', 1392])
        xps_table.append([67, 'Ho', '3', 'd', '5/2', 1351])
        xps_table.append([67, 'Ho', '4', 's', '0', 432.4])
        xps_table.append([67, 'Ho', '4', 'p', '1/2', 343.5])
        xps_table.append([67, 'Ho', '4', 'p', '3/2', 308.2])
        xps_table.append([67, 'Ho', '4', 'd', '3/2', 160])
        xps_table.append([67, 'Ho', '4', 'd', '5/2', 160])
        xps_table.append([67, 'Ho', '4', 'f', '5/2', 8.6])
        xps_table.append([67, 'Ho', '4', 'f', '7/2', 5.2])
        xps_table.append([67, 'Ho', '5', 's', '0', 49.3])
        xps_table.append([67, 'Ho', '5', 'p', '1/2', 30.8])
        xps_table.append([67, 'Ho', '5', 'p', '3/2', 24.1])

        xps_table.append([68, 'Er', '1', 's', '0', 57486])
        xps_table.append([68, 'Er', '2', 's', '0', 9751])
        xps_table.append([68, 'Er', '2', 'p', '1/2', 9264])
        xps_table.append([68, 'Er', '2', 'p', '3/2', 8358])
        xps_table.append([68, 'Er', '3', 's', '0', 2206])
        xps_table.append([68, 'Er', '3', 'p', '1/2', 2006])
        xps_table.append([68, 'Er', '3', 'p', '3/2', 1812])
        xps_table.append([68, 'Er', '3', 'd', '3/2', 1453])
        xps_table.append([68, 'Er', '3', 'd', '5/2', 1409])
        xps_table.append([68, 'Er', '4', 's', '0', 449.8])
        xps_table.append([68, 'Er', '4', 'p', '1/2', 366.2])
        xps_table.append([68, 'Er', '4', 'p', '3/2', 320.2])
        xps_table.append([68, 'Er', '4', 'd', '3/2', 167.6])
        xps_table.append([68, 'Er', '4', 'd', '5/2', 167.6])
        xps_table.append([68, 'Er', '4', 'f', '7/2', 4.7])
        xps_table.append([68, 'Er', '5', 's', '0', 50.6])
        xps_table.append([68, 'Er', '5', 'p', '1/2', 31.4])
        xps_table.append([68, 'Er', '5', 'p', '3/2', 24.7])

        xps_table.append([69, 'Tm', '1', 's', '0', 59390])
        xps_table.append([69, 'Tm', '2', 's', '0', 10116])
        xps_table.append([69, 'Tm', '2', 'p', '1/2', 9617])
        xps_table.append([69, 'Tm', '2', 'p', '3/2', 8648])
        xps_table.append([69, 'Tm', '3', 's', '0', 2307])
        xps_table.append([69, 'Tm', '3', 'p', '1/2', 2090])
        xps_table.append([69, 'Tm', '3', 'p', '3/2', 1885])
        xps_table.append([69, 'Tm', '3', 'd', '3/2', 1515])
        xps_table.append([69, 'Tm', '3', 'd', '5/2', 1468])
        xps_table.append([69, 'Tm', '4', 's', '0', 470.9])
        xps_table.append([69, 'Tm', '4', 'p', '1/2', 385.9])
        xps_table.append([69, 'Tm', '4', 'p', '3/2', 332.6])
        xps_table.append([69, 'Tm', '4', 'd', '3/2', 175.5])
        xps_table.append([69, 'Tm', '4', 'd', '5/2', 175.5])
        xps_table.append([69, 'Tm', '4', 'f', '7/2', 4.6])
        xps_table.append([69, 'Tm', '5', 's', '0', 54.7])
        xps_table.append([69, 'Tm', '5', 'p', '1/2', 31.8])
        xps_table.append([69, 'Tm', '5', 'p', '3/2', 25])

        xps_table.append([70, 'Yb', '1', 's', '0', 61332])
        xps_table.append([70, 'Yb', '2', 's', '0', 10486])
        xps_table.append([70, 'Yb', '2', 'p', '1/2', 9978])
        xps_table.append([70, 'Yb', '2', 'p', '3/2', 8944])
        xps_table.append([70, 'Yb', '3', 's', '0', 2398])
        xps_table.append([70, 'Yb', '3', 'p', '1/2', 2173])
        xps_table.append([70, 'Yb', '3', 'p', '3/2', 1950])
        xps_table.append([70, 'Yb', '3', 'd', '3/2', 1576])
        xps_table.append([70, 'Yb', '3', 'd', '5/2', 1528])
        xps_table.append([70, 'Yb', '4', 's', '0', 480.5])
        xps_table.append([70, 'Yb', '4', 'p', '1/2', 388.7])
        xps_table.append([70, 'Yb', '4', 'p', '3/2', 339.7])
        xps_table.append([70, 'Yb', '4', 'd', '3/2', 191.2])
        xps_table.append([70, 'Yb', '4', 'd', '5/2', 182.4])
        xps_table.append([70, 'Yb', '4', 'f', '5/2', 2.5])
        xps_table.append([70, 'Yb', '4', 'f', '7/2', 1.3])
        xps_table.append([70, 'Yb', '5', 's', '0', 52])
        xps_table.append([70, 'Yb', '5', 'p', '1/2', 30.3])
        xps_table.append([70, 'Yb', '5', 'p', '3/2', 24.1])

        xps_table.append([71, 'Lu', '1', 's', '0', 63314])
        xps_table.append([71, 'Lu', '2', 's', '0', 10870])
        xps_table.append([71, 'Lu', '2', 'p', '1/2', 10349])
        xps_table.append([71, 'Lu', '2', 'p', '3/2', 9244])
        xps_table.append([71, 'Lu', '3', 's', '0', 2491])
        xps_table.append([71, 'Lu', '3', 'p', '1/2', 2264])
        xps_table.append([71, 'Lu', '3', 'p', '3/2', 2024])
        xps_table.append([71, 'Lu', '3', 'd', '3/2', 1639])
        xps_table.append([71, 'Lu', '3', 'd', '5/2', 1589])
        xps_table.append([71, 'Lu', '4', 's', '0', 506.8])
        xps_table.append([71, 'Lu', '4', 'p', '1/2', 412.4])
        xps_table.append([71, 'Lu', '4', 'p', '3/2', 359.2])
        xps_table.append([71, 'Lu', '4', 'd', '3/2', 206.1])
        xps_table.append([71, 'Lu', '4', 'd', '5/2', 196.3])
        xps_table.append([71, 'Lu', '4', 'f', '5/2', 8.9])
        xps_table.append([71, 'Lu', '4', 'f', '7/2', 7.5])
        xps_table.append([71, 'Lu', '5', 's', '0', 57.3])
        xps_table.append([71, 'Lu', '5', 'p', '1/2', 33.6])
        xps_table.append([71, 'Lu', '5', 'p', '3/2', 26.7])

        xps_table.append([72, 'Hf', '1', 's', '0', 65351])
        xps_table.append([72, 'Hf', '2', 's', '0', 11271])
        xps_table.append([72, 'Hf', '2', 'p', '1/2', 10739])
        xps_table.append([72, 'Hf', '2', 'p', '3/2', 9561])
        xps_table.append([72, 'Hf', '3', 's', '0', 2601])
        xps_table.append([72, 'Hf', '3', 'p', '1/2', 2365])
        xps_table.append([72, 'Hf', '3', 'p', '3/2', 2107])
        xps_table.append([72, 'Hf', '3', 'd', '3/2', 1716])
        xps_table.append([72, 'Hf', '3', 'd', '5/2', 1662])
        xps_table.append([72, 'Hf', '4', 's', '0', 538])
        xps_table.append([72, 'Hf', '4', 'p', '1/2', 438.2])
        xps_table.append([72, 'Hf', '4', 'p', '3/2', 380.7])
        xps_table.append([72, 'Hf', '4', 'd', '3/2', 220])
        xps_table.append([72, 'Hf', '4', 'd', '5/2', 211.5])
        xps_table.append([72, 'Hf', '4', 'f', '5/2', 15.9])
        xps_table.append([72, 'Hf', '4', 'f', '7/2', 14.2])
        xps_table.append([72, 'Hf', '5', 's', '0', 64.2])
        xps_table.append([72, 'Hf', '5', 'p', '1/2', 38])
        xps_table.append([72, 'Hf', '5', 'p', '3/2', 29.9])

        xps_table.append([73, 'Ta', '1', 's', '0', 67416])
        xps_table.append([73, 'Ta', '2', 's', '0', 11682])
        xps_table.append([73, 'Ta', '2', 'p', '1/2', 11136])
        xps_table.append([73, 'Ta', '2', 'p', '3/2', 9881])
        xps_table.append([73, 'Ta', '3', 's', '0', 2708])
        xps_table.append([73, 'Ta', '3', 'p', '1/2', 2469])
        xps_table.append([73, 'Ta', '3', 'p', '3/2', 2194])
        xps_table.append([73, 'Ta', '3', 'd', '3/2', 1793])
        xps_table.append([73, 'Ta', '3', 'd', '5/2', 1735])
        xps_table.append([73, 'Ta', '4', 's', '0', 563.4])
        xps_table.append([73, 'Ta', '4', 'p', '1/2', 463.4])
        xps_table.append([73, 'Ta', '4', 'p', '3/2', 400.9])
        xps_table.append([73, 'Ta', '4', 'd', '3/2', 237.9])
        xps_table.append([73, 'Ta', '4', 'd', '5/2', 226.4])
        xps_table.append([73, 'Ta', '4', 'f', '5/2', 23.5])
        xps_table.append([73, 'Ta', '4', 'f', '7/2', 21.6])
        xps_table.append([73, 'Ta', '5', 's', '0', 69.7])
        xps_table.append([73, 'Ta', '5', 'p', '1/2', 42.2])
        xps_table.append([73, 'Ta', '5', 'p', '3/2', 32.7])

        xps_table.append([74, 'W', '1', 's', '0', 69525])
        xps_table.append([74, 'W', '2', 's', '0', 12100])
        xps_table.append([74, 'W', '2', 'p', '1/2', 11544])
        xps_table.append([74, 'W', '2', 'p', '3/2', 10207])
        xps_table.append([74, 'W', '3', 's', '0', 2820])
        xps_table.append([74, 'W', '3', 'p', '1/2', 2575])
        xps_table.append([74, 'W', '3', 'p', '3/2', 2281])
        xps_table.append([74, 'W', '3', 'd', '3/2', 1949])
        xps_table.append([74, 'W', '3', 'd', '5/2', 1809])
        xps_table.append([74, 'W', '4', 's', '0', 594.1])
        xps_table.append([74, 'W', '4', 'p', '1/2', 490.4])
        xps_table.append([74, 'W', '4', 'p', '3/2', 423.6])
        xps_table.append([74, 'W', '4', 'd', '3/2', 255.9])
        xps_table.append([74, 'W', '4', 'd', '5/2', 243.5])
        xps_table.append([74, 'W', '4', 'f', '5/2', 33.6])
        xps_table.append([74, 'W', '4', 'f', '7/2', 31.4])
        xps_table.append([74, 'W', '5', 's', '0', 75.6])
        xps_table.append([74, 'W', '5', 'p', '1/2', 45.3])
        xps_table.append([74, 'W', '5', 'p', '3/2', 36.8])

        xps_table.append([75, 'Re', '1', 's', '0', 71676])
        xps_table.append([75, 'Re', '2', 's', '0', 12527])
        xps_table.append([75, 'Re', '2', 'p', '1/2', 11959])
        xps_table.append([75, 'Re', '2', 'p', '3/2', 10535])
        xps_table.append([75, 'Re', '3', 's', '0', 2932])
        xps_table.append([75, 'Re', '3', 'p', '1/2', 2682])
        xps_table.append([75, 'Re', '3', 'p', '3/2', 2367])
        xps_table.append([75, 'Re', '3', 'd', '3/2', 1949])
        xps_table.append([75, 'Re', '3', 'd', '5/2', 1883])
        xps_table.append([75, 'Re', '4', 's', '0', 625.4])
        xps_table.append([75, 'Re', '4', 'p', '1/2', 518.7])
        xps_table.append([75, 'Re', '4', 'p', '3/2', 446.8])
        xps_table.append([75, 'Re', '4', 'd', '3/2', 273.9])
        xps_table.append([75, 'Re', '4', 'd', '5/2', 260.5])
        xps_table.append([75, 'Re', '4', 'f', '5/2', 42.9])
        xps_table.append([75, 'Re', '4', 'f', '7/2', 40.5])
        xps_table.append([75, 'Re', '5', 's', '0', 83])
        xps_table.append([75, 'Re', '5', 'p', '1/2', 45.6])
        xps_table.append([75, 'Re', '5', 'p', '3/2', 34.6])

        xps_table.append([76, 'Os', '1', 's', '0', 73871])
        xps_table.append([76, 'Os', '2', 's', '0', 12968])
        xps_table.append([76, 'Os', '2', 'p', '1/2', 12385])
        xps_table.append([76, 'Os', '2', 'p', '3/2', 10871])
        xps_table.append([76, 'Os', '3', 's', '0', 3049])
        xps_table.append([76, 'Os', '3', 'p', '1/2', 2792])
        xps_table.append([76, 'Os', '3', 'p', '3/2', 2457])
        xps_table.append([76, 'Os', '3', 'd', '3/2', 2031])
        xps_table.append([76, 'Os', '3', 'd', '5/2', 1960])
        xps_table.append([76, 'Os', '4', 's', '0', 658.2])
        xps_table.append([76, 'Os', '4', 'p', '1/2', 549.1])
        xps_table.append([76, 'Os', '4', 'p', '3/2', 470.7])
        xps_table.append([76, 'Os', '4', 'd', '3/2', 293.1])
        xps_table.append([76, 'Os', '4', 'd', '5/2', 278.5])
        xps_table.append([76, 'Os', '4', 'f', '5/2', 53.4])
        xps_table.append([76, 'Os', '4', 'f', '7/2', 50.7])
        xps_table.append([76, 'Os', '5', 's', '0', 84])
        xps_table.append([76, 'Os', '5', 'p', '1/2', 58])
        xps_table.append([76, 'Os', '5', 'p', '3/2', 44.5])

        xps_table.append([77, 'Ir', '1', 's', '0', 76111])
        xps_table.append([77, 'Ir', '2', 's', '0', 13419])
        xps_table.append([77, 'Ir', '2', 'p', '1/2', 12824])
        xps_table.append([77, 'Ir', '2', 'p', '3/2', 11215])
        xps_table.append([77, 'Ir', '3', 's', '0', 3174])
        xps_table.append([77, 'Ir', '3', 'p', '1/2', 2909])
        xps_table.append([77, 'Ir', '3', 'p', '3/2', 2551])
        xps_table.append([77, 'Ir', '3', 'd', '3/2', 2116])
        xps_table.append([77, 'Ir', '3', 'd', '5/2', 2040])
        xps_table.append([77, 'Ir', '4', 's', '0', 691.1])
        xps_table.append([77, 'Ir', '4', 'p', '1/2', 577.8])
        xps_table.append([77, 'Ir', '4', 'p', '3/2', 495.8])
        xps_table.append([77, 'Ir', '4', 'd', '3/2', 311.9])
        xps_table.append([77, 'Ir', '4', 'd', '5/2', 296.3])
        xps_table.append([77, 'Ir', '4', 'f', '5/2', 63.8])
        xps_table.append([77, 'Ir', '4', 'f', '7/2', 60.8])
        xps_table.append([77, 'Ir', '5', 's', '0', 95.2])
        xps_table.append([77, 'Ir', '5', 'p', '1/2', 63])
        xps_table.append([77, 'Ir', '5', 'p', '3/2', 48])

        xps_table.append([78, 'Pt', '1', 's', '0', 78395])
        xps_table.append([78, 'Pt', '2', 's', '0', 13880])
        xps_table.append([78, 'Pt', '2', 'p', '1/2', 13273])
        xps_table.append([78, 'Pt', '2', 'p', '3/2', 11564])
        xps_table.append([78, 'Pt', '3', 's', '0', 3296])
        xps_table.append([78, 'Pt', '3', 'p', '1/2', 3027])
        xps_table.append([78, 'Pt', '3', 'p', '3/2', 2645])
        xps_table.append([78, 'Pt', '3', 'd', '3/2', 2202])
        xps_table.append([78, 'Pt', '3', 'd', '5/2', 2122])
        xps_table.append([78, 'Pt', '4', 's', '0', 725.4])
        xps_table.append([78, 'Pt', '4', 'p', '1/2', 609.1])
        xps_table.append([78, 'Pt', '4', 'p', '3/2', 519.4])
        xps_table.append([78, 'Pt', '4', 'd', '3/2', 331.6])
        xps_table.append([78, 'Pt', '4', 'd', '5/2', 314.6])
        xps_table.append([78, 'Pt', '4', 'f', '5/2', 74.5])
        xps_table.append([78, 'Pt', '4', 'f', '7/2', 71.2])
        xps_table.append([78, 'Pt', '5', 's', '0', 101.7])
        xps_table.append([78, 'Pt', '5', 'p', '1/2', 65.3])
        xps_table.append([78, 'Pt', '5', 'p', '3/2', 51.7])

        xps_table.append([79, 'Au', '1', 's', '0', 80725])
        xps_table.append([79, 'Au', '2', 's', '0', 14353])
        xps_table.append([79, 'Au', '2', 'p', '1/2', 13734])
        xps_table.append([79, 'Au', '2', 'p', '3/2', 11919])
        xps_table.append([79, 'Au', '3', 's', '0', 3425])
        xps_table.append([79, 'Au', '3', 'p', '1/2', 3148])
        xps_table.append([79, 'Au', '3', 'p', '3/2', 2743])
        xps_table.append([79, 'Au', '3', 'd', '3/2', 2291])
        xps_table.append([79, 'Au', '3', 'd', '5/2', 2206])
        xps_table.append([79, 'Au', '4', 's', '0', 762.1])
        xps_table.append([79, 'Au', '4', 'p', '1/2', 642.7])
        xps_table.append([79, 'Au', '4', 'p', '3/2', 546.3])
        xps_table.append([79, 'Au', '4', 'd', '3/2', 353.2])
        xps_table.append([79, 'Au', '4', 'd', '5/2', 335.1])
        xps_table.append([79, 'Au', '4', 'f', '5/2', 87.6])
        xps_table.append([79, 'Au', '4', 'f', '7/2', 83.9])
        xps_table.append([79, 'Au', '5', 's', '0', 107.2])
        xps_table.append([79, 'Au', '5', 'p', '1/2', 74.2])
        xps_table.append([79, 'Au', '5', 'p', '3/2', 57.2])

        xps_table.append([80, 'Hg', '1', 's', '0', 83102])
        xps_table.append([80, 'Hg', '2', 's', '0', 14839])
        xps_table.append([80, 'Hg', '2', 'p', '1/2', 14209])
        xps_table.append([80, 'Hg', '2', 'p', '3/2', 12284])
        xps_table.append([80, 'Hg', '3', 's', '0', 3562])
        xps_table.append([80, 'Hg', '3', 'p', '1/2', 3279])
        xps_table.append([80, 'Hg', '3', 'p', '3/2', 2847])
        xps_table.append([80, 'Hg', '3', 'd', '3/2', 2385])
        xps_table.append([80, 'Hg', '3', 'd', '5/2', 2295])
        xps_table.append([80, 'Hg', '4', 's', '0', 802.2])
        xps_table.append([80, 'Hg', '4', 'p', '1/2', 680.2])
        xps_table.append([80, 'Hg', '4', 'p', '3/2', 576.6])
        xps_table.append([80, 'Hg', '4', 'd', '3/2', 378.2])
        xps_table.append([80, 'Hg', '4', 'd', '5/2', 358.8])
        xps_table.append([80, 'Hg', '4', 'f', '5/2', 104])
        xps_table.append([80, 'Hg', '4', 'f', '7/2', 99.9])
        xps_table.append([80, 'Hg', '5', 's', '0', 127])
        xps_table.append([80, 'Hg', '5', 'p', '1/2', 83.1])
        xps_table.append([80, 'Hg', '5', 'p', '3/2', 64.5])
        xps_table.append([80, 'Hg', '5', 'd', '3/2', 9.6])
        xps_table.append([80, 'Hg', '5', 'd', '5/2', 7.8])

        xps_table.append([81, 'Tl', '1', 's', '0', 85530])
        xps_table.append([81, 'Tl', '2', 's', '0', 15347])
        xps_table.append([81, 'Tl', '2', 'p', '1/2', 14698])
        xps_table.append([81, 'Tl', '2', 'p', '3/2', 12658])
        xps_table.append([81, 'Tl', '3', 's', '0', 3704])
        xps_table.append([81, 'Tl', '3', 'p', '1/2', 3416])
        xps_table.append([81, 'Tl', '3', 'p', '3/2', 2957])
        xps_table.append([81, 'Tl', '3', 'd', '3/2', 2485])
        xps_table.append([81, 'Tl', '3', 'd', '5/2', 2389])
        xps_table.append([81, 'Tl', '4', 's', '0', 846.2])
        xps_table.append([81, 'Tl', '4', 'p', '1/2', 720.5])
        xps_table.append([81, 'Tl', '4', 'p', '3/2', 609.5])
        xps_table.append([81, 'Tl', '4', 'd', '3/2', 405.7])
        xps_table.append([81, 'Tl', '4', 'd', '5/2', 385])
        xps_table.append([81, 'Tl', '4', 'f', '5/2', 122.2])
        xps_table.append([81, 'Tl', '4', 'f', '7/2', 117.8])
        xps_table.append([81, 'Tl', '5', 's', '0', 136])
        xps_table.append([81, 'Tl', '5', 'p', '1/2', 94.6])
        xps_table.append([81, 'Tl', '5', 'p', '3/2', 73.5])
        xps_table.append([81, 'Tl', '5', 'd', '3/2', 14.7])
        xps_table.append([81, 'Tl', '5', 'd', '5/2', 12.5])

        xps_table.append([82, 'Pb', '1', 's', '0', 88005])
        xps_table.append([82, 'Pb', '2', 's', '0', 15861])
        xps_table.append([82, 'Pb', '2', 'p', '1/2', 15200])
        xps_table.append([82, 'Pb', '2', 'p', '3/2', 13035])
        xps_table.append([82, 'Pb', '3', 's', '0', 3851])
        xps_table.append([82, 'Pb', '3', 'p', '1/2', 3554])
        xps_table.append([82, 'Pb', '3', 'p', '3/2', 3066])
        xps_table.append([82, 'Pb', '3', 'd', '3/2', 2586])
        xps_table.append([82, 'Pb', '3', 'd', '5/2', 2484])
        xps_table.append([82, 'Pb', '4', 's', '0', 891.8])
        xps_table.append([82, 'Pb', '4', 'p', '1/2', 761.9])
        xps_table.append([82, 'Pb', '4', 'p', '3/2', 643.5])
        xps_table.append([82, 'Pb', '4', 'd', '3/2', 434.3])
        xps_table.append([82, 'Pb', '4', 'd', '5/2', 412.2])
        xps_table.append([82, 'Pb', '4', 'f', '5/2', 141.7])
        xps_table.append([82, 'Pb', '4', 'f', '7/2', 136.9])
        xps_table.append([82, 'Pb', '5', 's', '0', 147])
        xps_table.append([82, 'Pb', '5', 'p', '1/2', 106.4])
        xps_table.append([82, 'Pb', '5', 'p', '3/2', 83.3])
        xps_table.append([82, 'Pb', '5', 'd', '3/2', 20.7])
        xps_table.append([82, 'Pb', '5', 'd', '5/2', 18.1])

        xps_table.append([83, 'Bi', '1', 's', '0', 90526])
        xps_table.append([83, 'Bi', '2', 's', '0', 16388])
        xps_table.append([83, 'Bi', '2', 'p', '1/2', 15711])
        xps_table.append([83, 'Bi', '2', 'p', '3/2', 13419])
        xps_table.append([83, 'Bi', '3', 's', '0', 3999])
        xps_table.append([83, 'Bi', '3', 'p', '1/2', 3696])
        xps_table.append([83, 'Bi', '3', 'p', '3/2', 3177])
        xps_table.append([83, 'Bi', '3', 'd', '3/2', 2688])
        xps_table.append([83, 'Bi', '3', 'd', '5/2', 2580])
        xps_table.append([83, 'Bi', '4', 's', '0', 939])
        xps_table.append([83, 'Bi', '4', 'p', '1/2', 805.2])
        xps_table.append([83, 'Bi', '4', 'p', '3/2', 678.8])
        xps_table.append([83, 'Bi', '4', 'd', '3/2', 464])
        xps_table.append([83, 'Bi', '4', 'd', '5/2', 440.1])
        xps_table.append([83, 'Bi', '4', 'f', '5/2', 162.3])
        xps_table.append([83, 'Bi', '4', 'f', '7/2', 157])
        xps_table.append([83, 'Bi', '5', 's', '0', 159.3])
        xps_table.append([83, 'Bi', '5', 'p', '1/2', 119])
        xps_table.append([83, 'Bi', '5', 'p', '3/2', 92.6])
        xps_table.append([83, 'Bi', '5', 'd', '3/2', 26.9])
        xps_table.append([83, 'Bi', '5', 'd', '5/2', 23.8])

        self.xps_table = list(xps_table)