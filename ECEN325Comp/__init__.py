import numpy as np
# Debug Message
def D(Message, Toggle=False):
  if Toggle:
    print(Message)

# Unit Manipulation
units = ["V", "A", 'Omega', "Ω", "F", '(A/V^2)']
Prefix = ["E", "P", "T", "G", "M", "k", "h", "da", "",
          "d", "c", "m", "u", "n", "p", "f", "a", 'SCI']
PrefixValue = [1E18, 1E15, 1E12, 1E9, 1E6, 1E3, 1E2, 1E1, 1,
                1E-1, 1E-2, 1E-3, 1E-6, 1E-9, 1E-12, 1E-15, 1E-18]

Standard_Prefixs = ["E", "P", "T", "G", "M", "k", "", "m", "u", "n", "p", "f", "a"]
Standard_PrefixValue = [1E18, 1E15, 1E12, 1E9, 1E6, 1E3, 1, 0.001, 1E-6, 1E-9, 1E-12, 1E-15, 1E-18]

def PrintUnits(v, v_unit, print_prefix=None, v_prefix=None, RoundingDigits=2,debug=False):
  # Major Errors:
  if v == None or v_unit == None:
    return "[Undefined]"
  elif v > PrefixValue[0]:
    return f"∞"
  elif v == 0:
    return f"0 {v_unit}"
  # Auto-Find print-prefix if unknown
  if print_prefix == 'SCI' and v_prefix == None:
    D(f"v")
    for i in range(len(Standard_Prefixs)):
      D(f"v = {v}, i = {Standard_Prefixs[i]}, v > i\n   index = {i}", debug)
      if abs(v) > Standard_PrefixValue[i]:
        D(f"v = {v}, i = {Standard_PrefixValue[i]}, v > i", debug)
        print_prefix = Standard_Prefixs[i]
        break
    if print_prefix == None or print_prefix == 'SCI':
      raise Exception(f"Invalid Values:\n V = {v}\nv_unit={v_unit}\nPrintPrefix = {print_prefix}")
  

  # Redefine the Units:
  if v_unit.lower() == "Omega".lower():
    unit = "Ω"
  else:
    unit = v_unit
  # Print Units if no conversion
  if v_prefix == None and print_prefix == None:
    return f"{round(v, RoundingDigits)} {unit}"


  # --- --- ---
  # Find Basic Stuff
  if v_prefix == None:
    v_prefix = ""
  if unit in units and v_prefix in Prefix and print_prefix in Prefix:
    v_index = Prefix.index(v_prefix)
    print_index = Prefix.index(print_prefix)
    vi = v * PrefixValue[v_index]
    vp = vi / PrefixValue[print_index]
    D(f"PrefixValue[v_index] = {PrefixValue[v_index]} [-] vi = {vi}, vp = {vp}", debug)
    return f"{round(vp, RoundingDigits)} {print_prefix}{unit}"


  # --- --- ---
  # Error Detection
  if not unit in units:
    raise Exception(f"Invalid Unit {unit}")
  elif not v_prefix in Prefix:
    raise Exception(f"Invalid Value Prefix: {v_prefix}")
  elif not print_prefix in Prefix:
    raise Exception(f"Invalid print prefix {print_prefix}")


# Operations:
def parallel(ListOfValues):
  InverseSum = 0
  for i in ListOfValues:
    InverseSum += (i ** -1)
  return (InverseSum ** -1)

# MOSFET Simulator
# @title MOS Class
# @title MOS Class
class MOS:
  """
  This class is meant to simulate a MOSFET transister.

  NAME(str): This is the name of the string and will appear be used when the
              circuit attributes are printed
  TYPE(str): This is the type of MOSFET transister. Must be either 'PMOS' or
              'NMOS'
  NHigh(str): This is the terminal with the highest voltage level. Necessary for
              calculating Voltage Drop and Voltage at terminals.
  IN(char): This is the input resistance for the MOSFET. Must be 's' (source),
             'g' (gate), or 'd' (drain).
  OUT(char): This is the output resistance for the MOSFET. Must be 's' (source),
             'g' (gate), or 'd' (drain).
  NOTE: EITHER MOSBETA or Width, Length, and KPRIME MUST be provided or
              circuit analysis is impossible.
  Width(float): This is the physical width of the MOSFET. Used in finding
                  MOSBETA and the current when in a current mirror.
  Length(float): This is the physical length of the MOSFET. Used in finding
                  MOSBETA and the current when in a current mirror.
  KPRIME(float): This is the Transconductance Parameter of the MOSFET. It
                  varies based on the physical makeup and fabrication process
                  of the MOSFET. Used for calculating MOSBETA.
  MOSBETA(float): Value is determined by KPRIME * (Width / Length). Used for
                  calculating MOSFET current.
  Vt(float): The Threshold Voltage of a MOSFET. Used for determining the
              region of the device.
  MirrorBase(int): This value indicates whether the MOSFET is current mirrored.
                    Possible States:
                    0 (default). Current is completely independant. 
                    1. Current Dependant. Current is dependant on the current 
                    through the base MOSFET. This MOSFET is in series with 
                    either a Current Mirror or Current Mirror Base. 
                    2. Current Mirror: The current through this MOSFET is the
                    same as through the Mirror Base, and it shares a node with 
                    the node VG. 
                    3. This MOSFET is a Mirror Base. The Drain and Gate
                    teminals are connect and share the same voltage. VDG is set
                    to 0. 
  MirrorList([MOS]): Only used if MirrorBase == 3. This is a list of all the 
                      MOSFETS that are updated whenever the voltage on the 
                      Mirror Base is updated. 

  """
  def __init__(self, NAME, TYPE, NHigh=None, IN=None, OUT=None, 
               VS=None, VG=None, VD=None,
               Vt = None, 
               Width=1, Length=1, 
               KPRIME = None, MOSBETA = None, 
               MirrorBase=0, dependant_list=None, mirror_list=None, 
               debug=False):
    # Define values that MUST be correct:
    self._NAME = NAME
    self._DEBUG = debug
    self._TYPE = TYPE if TYPE in ['PMOS', 'NMOS'] else print(f"MOS {self._NAME} Type is invalid. MUST BE (PMOS/NMOS). Please Update")
    if IN != None and OUT != None:
      self._ri = IN if IN in ['s', 'g', 'd'] else print(f"MOS {self._NAME} Input Resistance is invalid. MUST BE (s/g/d). Please Update")
      self._ro = OUT if OUT in ['s', 'g', 'd'] else print(f"MOS {self._NAME} Output Resistance is invalid. MUST BE (s/g/d). Please Update")
    else:
      self._ri = IN
      self._ro = OUT
    # Define Fundamental Elements

    # Component Values:
    self._Vt = Vt
    self._KPRIME = KPRIME
    self._W = Width
    self._L = Length
    self._MOSBETA = MOSBETA
    self.CALC_MOSBETA()
    # Check For Activation:
    self._SAT = None
    self._CUTOFF = None
    # --- --- ---
    # Define Voltage At Nodes:
    self._NHigh = NHigh
    self._VD = VD
    self._VG = VG
    self._VS = VS
    if MirrorBase == 3:
      self._VLIST = [self._VD, self._VG]
      self._VG = self._VLIST[0]
      self._VD = self._VLIST[0]
    # Define Voltage Differences:
    self._VDS = None
    self._VDG = None
    self._VGS = None
    self._Vov = None
    # Define the current:
    self._MIRRORBASE = MirrorBase
    self._CURRENT_DEPENDANTS = dependant_list
    self._CURRENT_MIRRORS = mirror_list
    self._ID = None
    # Define values for solving the quadratic to find VGS and ID:
    self.VI_CALCTYPE = None
    self._VGSp = None
    self._VGSm = None
    # --- --- ---
    # Define AC Values:
    # Define Small Signal Model Components
    self._gm = None

    # Define External Resistances
    self._ZS = None
    self._ZG = None
    self._ZD = None

    # Define Various Input Resistances
    self._Z_source = None
    self._Z_gate = float('inf')
    self._Z_drain = float('inf')
  

  def NODESET(self, VD=None, VG=None, VS=None):
    if VD != None:
      self._VD = VD
    if VG != None:
      self._VG = VG
    if self._VS != None:
      self._VS = VS
  # --- --- --- Print MOS Metrics --- --- ---

  def PrintAttributes(self):
    '''
    This function prints all of the attributes of this MOSFET in a disorganized
    list. Very Useful. 
    '''
    a = '   '
    print(f" --- --- --- \n"
          f"Metrics for MOS {self._NAME}\n"
          f"{a}Highest Voltage Node = {self._NHigh}\n"
          f"{a}Type = {self._TYPE} [-] Vtn = {PrintUnits(self._Vt, 'V', 'SCI')} [-] β = {PrintUnits(self._MOSBETA, '(A/V^2)', 'SCI', debug=False)}\n"
          f"DC Solutions:\n"
          f"{a}VD = {PrintUnits(self._VD, 'V', 'SCI')} [-] VS = {PrintUnits(self._VS, 'V', 'SCI')} [-] VG = {PrintUnits(self._VG, 'V', 'SCI')}\n"
          f"{a}VDS = {PrintUnits(self._VDS, 'V', 'SCI')} [-] VGS = {PrintUnits(self._VGS, 'V', 'SCI')} [-] VDG = {PrintUnits(self._VDG, 'V', 'SCI')}\n"
          f"{a}Vov = {PrintUnits(self._Vov, 'V', 'SCI')} [-] ID = {PrintUnits(self._ID, 'A', 'SCI')}"
          f" --- --- --- \n"
          f"AC Solutions:\n"
          f"{a}ro = {PrintUnits(self.ro, 'Omega', 'SCI')} [-] ri = {PrintUnits(self.ri, 'Omega', 'SCI')}\n"
          f"{a}gm = {PrintUnits(self._gm, 'Omega', 'SCI')}\n"
          f"{a}ZS = {PrintUnits(self._ZS, 'Omega', 'SCI', debug=self._DEBUG)} [-] ZG = {PrintUnits(self._ZG, 'Omega', 'SCI', debug=self._DEBUG)} [-] ZD = {PrintUnits(self._ZD, 'Omega', 'SCI', debug=self._DEBUG)}\n"
          f"{a}ZSource = {PrintUnits(self._Z_source, 'Omega', 'SCI', debug=self._DEBUG)} [-] ZGate = {PrintUnits(self._Z_gate, 'Omega', 'SCI', debug=self._DEBUG)} [-] ZDrain = {PrintUnits(self._Z_drain, 'Omega', 'SCI', debug=self._DEBUG)}")
    if self.VI_CALCTYPE:
      print(f"Parallel Voltage Method was used.\n"
            f"{a}VGSp = {PrintUnits(self._VGSp, 'V', 'SCI')} [-] VGSm = {PrintUnits(self._VGSm, 'V', 'SCI')}\n"
            f"{a}VGS = {PrintUnits(self._VGS, 'V', 'SCI')}")
    self.Active()

  def Active(self):
    D("'Active' Function Called.", self._DEBUG)
    '''
    This functions prints what active region the MOSFET is in with the current
    values.
    '''
    a = '   '
    if (self._SAT == None or self._CUTOFF == None):
      self.CheckRegion()    
    print(f"MOSFET {self._NAME} Region:")
    if self._VGS != None and self._Vt != None and self._CUTOFF:
      print(f"{a}VGS = {PrintUnits(self._VGS, 'V', 'SCI')}, Vtn = {PrintUnits(self._Vt, 'V', 'SCI')}\n"
            f"{a}VGS < Vtn\n"
            f"{a}MOS is in Cutoff Region.")
    elif self._VDS != None and self._Vov != None and self._SAT:
      print(f"{a}VDS = {PrintUnits(self._VDS, 'V', 'SCI')}, Vov = {PrintUnits(self._Vov, 'V', 'SCI')}\n"
            f"{a}{PrintUnits(self._VDS, 'V', 'SCI')} > {PrintUnits(self._Vov, 'V', 'SCI')}\n"
            f"{a}MOS is in active region.")
    elif self._VDS != None and self._Vov != None and not self._SAT:
      print(f"{a}VDS = {PrintUnits(self._VDS, 'V', 'SCI')}, Vov = {PrintUnits(self._Vov, 'V', 'SCI')}\n"
            f"{a}{PrintUnits(self._VDS, 'V', 'SCI')} < {PrintUnits(self._Vov, 'V', 'SCI')}\n"
            f"{a}MOS is in Triode Region.")
    else:
      print(f"VDS = {PrintUnits(self._VDS, 'V', 'SCI')}, Vov = {PrintUnits(self._Vov, 'V', 'SCI')}\n"
            f"MOS region is unknown.")
    print("\n")



  # --- --- --- Check whether is Active region --- --- ---

  def CheckRegion(self):
    '''
    This function checks whether the MOSFET is in the CUTOFF or SAT region. 
    If VSG or Vt are undefined, then CUTOFF cannot be calculated. 
    If VDS or Vov are undefined, then SATURATION cannot be calculated.
    '''
    D("CheckRegion Function Called", self._DEBUG)
    if self._CUTOFF == None and self._VGS != None and self._Vt != None:
      self._CUTOFF = self._VGS < self._Vt
      D(f"[CheckRegion Function set variable self._CUTOFF = {self._CUTOFF}]", self._DEBUG)
    if self._SAT == None and self._VDS != None and self._Vov != None:
      self._SAT = self._VDS > self._Vov
      D(f"[CheckRegion Function set variable self._SAT = {self._SAT}]", self._DEBUG)


  # --- --- --- Calculate Values Basic Values --- --- ---

  def CALC_CURRENTMIRROR(self):
    '''
    This function changes the values of all MOSFETS in the _CURRENT_MIRROR and 
    _CURRENT_DEPENDANT lists.
    '''
    if self._MIRRORBASE == 3:
      for MOS in self._CURRENT_MIRRORS:
        if self._ID != None:
          MOS.ID = self._ID
        if self._VG != None:
          MOS.VG = self._VG
      
      for MOS in self._CURRENT_DEPENDANTS:
        if self._ID != None:
          MOS.ID = self._ID


  # --- --- --- Calculate Changed Values for MOSBETA --- --- ---

  def CALC_MOSBETA(self):
    D("CALC_MOSBETA Function Called", self._DEBUG)
    '''
    Calculates β of the MOSFET when k', W, and L are defined
    within the class. 
    
    Will do nothing if β (_MOSBETA) is already defined, 
    or if any of the necessary variables are undefined.
    '''
    '_MOSBETA depends on L, W, and kprime'
    if self._MOSBETA == None and self._KPRIME != None and self._W != None and self._L != None:
      self._MOSBETA = self._KPRIME * (self._W / self._L)
      D("CALC_MOSBETA set variable MOSBETA", self._DEBUG)
    if self._KPRIME == None and self._MOSBETA != None and self._W != None and self._L != None:
      self._KPRIME = self._MOSBETA * (self._L / self._W)


  def CALC_Vov(self):
    
    '''
    1. if VGS and Vt are known, then it calculates Vov (Overdrive Voltage).
    2. If Vov or ID are known, then it calculates gm (transconductance)
    3. If Vov is known, then it runs the function CheckRegion() which checks
    what region the MOSFET is in.
    '''
    D("CALC_Vov Function Called", self._DEBUG)
    # - - - - - 
    if self._Vov == None and self._VGS != None and self._Vt != None:
      self._Vov = abs(self._VGS) - abs(self._Vt)
      D("CALC_Vov set variable Vov", self._DEBUG)
    if self._Vov != None:
      self.CheckRegion()
    

  def CALC_VDS(self):
    '''
    Calculates VDS if VD and VS are known.
    '''
    D("CALC_VDS function called.", self._DEBUG)
    if self._VD != None and self._VS != None and self._VDS == None:
      D("CALC_VDS set variable VDS", self._DEBUG)
      self._VDS = abs(self._VD - self._VS)
      self.CheckRegion()


  def CALC_VGS(self):
    '''
    Calculates VGS if VG and VS are known.
    '''
    D("CALC_VGS function called.", self._DEBUG)
    if self._VG != None and self._VS != None and self._VGS == None:
      D("CALC_VGS set variable VGS", self._DEBUG)
      self._VGS = abs(self._VG - self._VS)
      self.CALC_CURRENTMIRROR()


  def CALC_VDG(self):
    '''
    Calculates VDG if VD and VG are known.
    '''
    D("CALC_VDG function called.", self._DEBUG)
    if self._VD != None and self._VG != None and self._VDG == None:
      D("CALC_VDG set variable VDG", self._DEBUG)
      self._VDG = abs(self._VD - self._VG)
      self.CALC_CURRENTMIRROR()
  

  def CALC_GM(self):
    '''
    This function calculates gm (transconductance) if β (MOSBETA) is known and 
    either Vov (overdrive voltage) or ID (current) are known
    '''
    D("CALC_GM function called.", self._DEBUG)
    if self._Vov != None and self._MOSBETA != None:
      self._gm = self._MOSBETA * self._Vov
      D("CALC_GM set variable gm (1)", self._DEBUG)
    elif self._ID != None and self._MOSBETA != None:
      self._gm = (2 * self._MOSBETA * self._ID) ** (1/2)
      D("CALC_GM set variable gm (2)", self._DEBUG)
    if self._gm != None:
      self._Z_source = 1 / self._gm
  

  def CALC_NHIGH(self):
    '''
    If NHigh isn't given in the class definition, this function attempts to 
    find NHigh. The voltage at two terminals is required otherwise nothing will
    happen.
    '''
    D("CALC_NHIGH function called.", self._DEBUG)
    if self._NHigh != None:
      D("NHigh variable is already known. Returning.", self._DEBUG)
      return
    # Comparison of VS and VG
    if self._VS != None and self._VG != None:
      if self._VG > self._VS:
        # If the Gate is above the source, then the drain is above the source.
        self._NHigh = 'd'
        D(f"CALC_NHIGH set variable NHigh = {self._NHigh}.", self._DEBUG)
        return
      if self._VG < self._VS:
        # If the Gate is below the source, then the drain is below the source
        self._NHigh = 's'
        D(f"CALC_NHIGH set variable NHigh = {self._NHigh}.", self._DEBUG)
        return
    

    # Comparison of VS and VD
    if self._VS != None and self._VD != None:
      if self._VD > self._VS:
        # If the Drain is above the source. 
        self._NHigh = 'd'
        D(f"CALC_NHIGH set variable NHigh = {self._NHigh}.", self._DEBUG)
        return
      if self._VD < self._VS:
        # If the Source is above the drain. 
        self._NHigh = 's'
        D(f"CALC_NHIGH set variable NHigh = {self._NHigh}.", self._DEBUG)
        return
    
    # Comparison of VD and VG
    if self._VG != None and self._VD != None:
      if self._VG > self._VD:
        # If the Drain is below the gate, then the source is above the drain.
        self._NHigh = 's'
        D(f"CALC_NHIGH set variable NHigh = {self._NHigh}.", self._DEBUG)
        return
      if self._VG < self._VD:
        # If the Drain is above the gate, then the source is below the Drain.
        self._NHigh = 'd'
        D(f"CALC_NHIGH set variable NHigh = {self._NHigh}.", self._DEBUG)
        return

  # --- --- --- Calculate Voltage Node Voltage From Voltage Difference --- --- ---

  def CALC_VS_HIGH(self):
    '''
    Calculates various voltage relationships from known values. Only runs if
    NHigh = 's' (Source).
    '''
    D("CALC_VS_HIGH function called.", self._DEBUG)
    # Check for Change in VGS
    D("CALC_VS_HIGH Checking for Change in VGS", self._DEBUG)
    if self._VGS != None and self._VS != None and self._VG == None:
      self._VG = self._VS - self._VGS
    if self._VGS != None and self._VS == None and self._VG != None:
      self._VS = self._VG + self._VGS
    
    # Check for Change in VDG
    D("CALC_VS_HIGH Checking for Change in VDG", self._DEBUG)
    if self._VDG != None and self._VG != None and self._VD == None:
      self._VD = self._VG - self._VDG
    if self._VDG != None and self._VG == None and self._VD != None:
      self._VG = self._VD + self._VDG
    
    # Check for change in VDS
    D("CALC_VS_HIGH Checking for Change in VDS", self._DEBUG)
    if self._VDS != None and self._VS == None and self._VD != None:
      self._VD = self._VS - self._VDS
    if self._VDS != None and self._VS != None and self._VD == None:
      self._VS = self._VD + self._VDS
  # ---


  def CALC_VD_HIGH(self):
    '''
    Calculates various voltage relationships from known values. Only runs if
    NHigh = 'd' (Drain).
    '''
    D("CALC_VD_HIGH function called.", self._DEBUG)
    # Check for Change in VGS
    D("CALC_VD_HIGH Checking for Change in VGS", self._DEBUG)
    if self._VGS != None and self._VS != None and self._VG == None:
      self._VG = self._VS + self._VGS
    if self._VGS != None and self._VS == None and self._VG != None:
      self._VS = self._VG - self._VGS
    
    # Check for Change in VDG
    D("CALC_VD_HIGH Checking for Change in VDG", self._DEBUG)
    if self._VDG != None and self._VG != None and self._VD == None:
      self._VD = self._VG + self._VDG
    if self._VDG != None and self._VG == None and self._VD != None:
      self._VG = self._VD - self._VDG
    
    # Check for change in VDS
    D("CALC_VD_HIGH Checking for Change in VDS", self._DEBUG)
    if self._VDS != None and self._VS == None and self._VD != None:
      self._VD = self._VS + self._VDS
    if self._VDS != None and self._VS != None and self._VD == None:
      self._VS = self._VD - self._VDS
  # ---


  def CALC_VNODE(self):
    '''
    Calculates various voltage relationships from known values. Only runs if 
    NHigh is known. 
    '''
    D("CALC_VNODE function called.", self._DEBUG)
    if self._NHigh == None:
      D("CALC_VNODE called CALC_NHIGH", self._DEBUG)
      self.CALC_NHIGH()
      if self._NHigh == None:
        D("NHigh Still Unknown. Exiting self CALC_VNODE. ", self._DEBUG)
        return
    if self._NHigh == 's':
      D("", self._DEBUG)
      self.CALC_VS_HIGH()
      self.CALC_VS_HIGH()
    elif self._NHigh == 'd':
      D("", self._DEBUG)
      self.CALC_VD_HIGH()
      self.CALC_VD_HIGH()
    self.CALC_Vov()  
    self.CALC_CURRENTMIRROR()


  # --- --- --- Unique Parameter Calculation Methods --- --- ---

  def SeriesMirror(self, M, VH=5, VL=0, Resistance=0, r=2):
    '''
    DO NOT USE WITHOUT UNDERSTANDING!!!
    This function calculates the VSG across two seperate Current Mirror Bases.
    '''
    D(f"Series Mirror Function Called: ")
    ΔV = abs(VH - VL)
    D(f"Resistace = {Resistance}", self._DEBUG or M._DEBUG)
    M1 = self
    M2 = M
    if self.VS > M.VS:
      Y1 = (M1._MOSBETA / M2._MOSBETA) ** (1/2)
      Z1 = Y1 * abs(M1._Vt) - abs(M2._Vt)
      Y2 = 1 / Y1
      Z2 = Z1 / Y1
      # --- --- ---
      # Setup and solve quadratic for M1 == self
      A1 = M1._MOSBETA * Resistance * (1/2)
      B1 = (-1 * Resistance * M1._MOSBETA * abs(M1._Vt)) + (Y1 + 1)
      C1 = (M1._Vt ** 2) * (M1._MOSBETA) * (Resistance) * (1/2) - ΔV + Z1
      D(f"A1 = {A1}, B1 = {B1}, C1 = {C1}", self._DEBUG)
      self._VGSp = ( (-1 * B1) + (((B1**2) - 4*A1*C1) ** (1/2)) ) / (2 * A1)
      self._VGSm = ( (-1 * B1) - (((B1**2) - 4*A1*C1)  ** (1/2)) ) / (2 * A1)
      D(f"{M1._NAME}.VGSp = {M1._VGSp} [-] {M1._NAME}.VGSm = {M1._VGSm}", M1._DEBUG)

      # Setup and solve quadratic for M2
      A2 = M2._MOSBETA * Resistance * (1/2)
      B2 = (-1 * Resistance * M2._MOSBETA * abs(M2._Vt)) + (Y2 + 1)
      C2 = (M2._Vt ** 2) * (M2._MOSBETA) * (Resistance) * (1/2) - ΔV + Z2
      D(f"A2 = {A2}, B2 = {B2}, C2 = {C2}", M2._DEBUG)
      M2._VGSp = ( (-1 * B2) + (((B2**2) - 4*A2*C2) ** (1/2)) ) / (2 * A2)
      M2._VGSm = ( (-1 * B2) - (((B2**2) - 4*A2*C2) ** (1/2)) ) / (2 * A2)
      D(f"{M2._NAME}.VGSp = {M2._VGSp} [-] {M2._NAME}.VGSm = {M2._VGSm}", M2._DEBUG)

      # Check for which value of VGS is value for M1 == self:
      if self._VGSp < abs(self._Vt) and self._VGSm >= self._Vt:
        self._VGS = self._VGSm
      elif self._VGSp >= abs(self._Vt) and self._VGSm < abs(self._Vt):
        self._VGS = self._VGSp
      elif self._VGSp < abs(self._Vt) and self._VGSm < abs(self._Vt):
        raise Exception(f"VGSp = {self._VGSp} [-] VGSn = {self._VGSm}\nINVALID {self._NAME}.VGS value. CODE: 1")
      else:
        raise Exception(f"VGSp = {self._VGSp} [-] VGSn = {self._VGSm}\nINVALID {self._NAME}.VGS value. CODE: 2")
      
      
      # Check Value for M2:
      if M2._VGSp < abs(M2._Vt) and M2._VGSm >= M2._Vt:
        M2._VGS = M2._VGSm
      elif M2._VGSp >= abs(M2._Vt) and M2._VGSm < abs(M2._Vt):
        M2._VGS = M2._VGSp
      elif M2._VGSp < abs(M2._Vt) and M2._VGSm < abs(M2._Vt):
        raise Exception(f"VGSp = {M2._VGSp} [-] VGSn = {M2._VGSm}\nINVALID {M2._NAME}.VGS value. CODE: 1")
      else:
        raise Exception(f"VGSp = {M2._VGSp} [-] VGSn = {M2._VGSm}\nINVALID {M2._NAME}.VGS value. CODE: 2")

      
      # Calculate ID for both. 
      M1.ID = M2.ID = (ΔV - M1.VGS - M2.VGS) * (1 / Resistance)
      self.CALC_GM()
      self.CALC_VNODE()
      self.CALC_Vov()
      M2.CALC_GM()
      M2.CALC_VNODE()
      M2.CALC_Vov()



      



  def ParallelVoltageMethod(self, VP, R, r=4, CCI=None):
    '''
    DO NOT USE WITHOUT UNDERSTANDING!!!
    This function calculates the voltage VGS and current ID when given a
    parallel voltage VP, and a resistance R. If the voltage going through the
    resistor R is from more than one source, CCI is current from other sources
    traveling through resistor R.
    '''
    D("", self._DEBUG)
    # Check for CCI change
    if CCI != None:
      VP = VP - R * CCI
      # print(R6.V - (R4.Z * M1._ID))
    # First, we check whether Vov is known. Must know VGS and Vt
    if self._ID == None and self._Vov != None and self._MOSBETA != None:
      D("", self._DEBUG)
      self.VI_CALCTYPE = False
      self._ID = (self._MOSBETA * (self._Vov ** 2))/(2)
      self.CALC_GM()


    # Now check
    elif self._ID == None and self._VGS == None and self._Vt != None:
      self.VI_CALCTYPE = True
      D("", self._DEBUG)


      # Find the parts of the equation
      A = self._MOSBETA / 2
      B = ((1/R) - (self._MOSBETA * abs(self._Vt)))
      C = (self._MOSBETA / 2) * (abs(self._Vt) ** 2) - (abs(VP) /R)
      if self._DEBUG:
        print(f"A = {A}, B = {B}, C = {C}")


      # Now use the Quadratic Formula
      self._VGSp = round( ( (-1 * B) + (((B**2) - 4*A*C)) ** (1/2)) / (2 * A), r)
      self._VGSm = round( ( (-1 * B) - (((B**2) - 4*A*C)) ** (1/2)) / (2 * A), r)


      # Now to determine which value to use.
      if self._VGSp < abs(self._Vt) and self._VGSm >= self._Vt:
        self._VGS = self._VGSm
      elif self._VGSp >= abs(self._Vt) and self._VGSm < abs(self._Vt):
        self._VGS = self._VGSp
      elif self._VGSp < abs(self._Vt) and self._VGSm < abs(self._Vt):
        raise Exception(f"VGSp = {self._VGSp} [-] VGSn = {self._VGSm}\nINVALID {self._NAME}.VGS value. CODE: 1")
      else:
        raise Exception(f"VGSp = {self._VGSp} [-] VGSn = {self._VGSm}\nINVALID {self._NAME}.VGS value. CODE: 2")
      

      # Now reverse calculate the current ID:
      if CCI == None:
        _CCI = 0
      else:
        _CCI = CCI
      IR1 = (VP - self._VGS) / R
      self._ID = IR1
      self.CALC_GM()
      
      
      # Finally, call VNODE to update VG or VS as necessary.
      self.CALC_VNODE()
      self.CALC_Vov()
      

  # --- --- --- Define all the properties and setters --- --- ---
  # --- Define DC Voltage Functions ---
  

  @property
  def VD(self):
    '''
    This function retrieves the value of the voltage at the Drain terminal of 
    the MOSFET. 
    '''
    D("", self._DEBUG)
    if self._MIRRORBASE == 3:
      return self._VG
    return self._VD

  @VD.setter
  def VD(self, v):
    '''
    This function sets the value of the voltage at the Drain terminal of the
    MOSFET. 
    
    '''
    D("", self._DEBUG)
    self._VD = v
    if self._MIRRORBASE == 3:
      self._VG = v
    # Define Related Values
    self.CALC_VDS()
    self.CALC_VDG()
    self.CALC_VNODE()
  @property
  def VG(self):
    '''
    Returns the voltage at the Gate terminal of the MOSFET
    '''
    D("", self._DEBUG)
    return self._VG
  @VD.setter
  def VG(self, v):
    '''
    Sets the voltage on the Gate terminal of the MOSFET. 
    Also calculates related values.
    '''
    D("", self._DEBUG)
    self._VG = v
    if self._MIRRORBASE == 3:
      self._VD = v
      self.CALC_VDS()
      self.CALC_VDG()
    # Calculate any related values
    self.CALC_VGS()
    self.CALC_VDG()
    self.CALC_VNODE()

  @property
  def VS(self):
    '''
    Returns the voltage at the Source terminal of the MOSFET
    '''
    D("", self._DEBUG)
    return self._VS
  @VS.setter
  def VS(self, v):
    '''
    Sets the voltage on eht Source Terminal. Also calculates related values.
    '''
    D("", self._DEBUG)
    self._VS = v
    self.CALC_VDS()
    self.CALC_VGS()
    self.CALC_VNODE()

  # --- Define Voltage Difference Functions ---

  @property
  def VDS(self):
    '''
    Returns the voltage difference between the Drain and Source nodes. 
    '''
    D("", self._DEBUG)
    return self._VDS

  @property
  def VGS(self):
    '''
    Returns the voltage difference between the Gate and Source nodes.
    '''
    D("", self._DEBUG)
    return self._VGS

  @property
  def VDG(self):
    '''
    Returns the voltage difference between the Gate and Drain nodes.
    '''
    D("", self._DEBUG)
    return self._VDG

  # --- Define Current Functions ---

  @property
  def ID(self):
    '''
    Returns the current through the MOSFET. 
    '''
    D("Property ID called.", self._DEBUG)
    return self._ID
  @ID.setter
  def ID(self, v):
    '''
    Sets the current through the MOSFET. 
    Calculates Vov (Overdrive Voltage) when β is known. 
    Calculates gm (Tranconductance). 
    '''
    D("ID Setter called.", self._DEBUG)
    self._ID = v
    self.CALC_Vov()
    self.CALC_GM()
    self.CALC_CURRENTMIRROR()
  @property
  def I(self):
    return self.ID
  @I.setter
  def I(self, v):
    self.ID = v
  
  @property
  def dlist(self):
    return self._CURRENT_DEPENDANTS
  @dlist.setter
  def dlist(self, v):
    self._CURRENT_DEPENDANTS = v
  
  @property
  def mlist(self):
    return self._CURRENT_MIRRORS
  @mlist.setter
  def mlist(self,v):
    self._CURRENT_MIRRORS = v

  # --- --- --- Define all the AC Values --- --- ---
  # Define Interal AC Resistances:
  @property
  def ZSource(self):
    '''
    Returns the input resistance through the source terminal.
    '''
    D("Property ZSource called.", self._DEBUG)
    return self._Z_source
  @property
  def ZGate(self):
    '''
    Returns the input resistance through the gate terminal.
    '''
    D("Property ZGate called.", self._DEBUG)
    return self._Z_gate
  @property
  def ZDrain(self):
    '''
    Returns the input resistance through the drain terminal.
    '''
    D("Property ZDrain called.", self._DEBUG)
    return self._Z_drain
    
  # Define External AC Resistances
  @property
  def Z(self):
    '''

    '''
    D("", self._DEBUG)
    return self._Z
    
  @Z.setter
  def ZS(self, v):
    '''

    '''
    D("", self._DEBUG)
    self._Z = v
  

  @property
  def ZS(self):

    D("", self._DEBUG)
    return self._ZS
  @ZS.setter
  def ZS(self, v):

    D("", self._DEBUG)
    self._ZS = v
  

  @property
  def ZG(self):

    D("", self._DEBUG)
    return self._ZG
  @Z.setter
  def ZG(self, v):

    D("", self._DEBUG)
    self._ZG = v
  

  @property
  def ZD(self):

    D("", self._DEBUG)
    return self._ZD
  @Z.setter
  def ZD(self, v):

    D("", self._DEBUG)
    self._ZD = v

  # Properties of "ro" Input resistance
  @property
  def ri(self):
    '''
    Returns the input resistance along the defined path. 
    '''
    D("", self._DEBUG)
    if self._ri == 'd':
      return self._Z_drain
    elif self._ri == 'g':
      return self._Z_gate
    elif self._ri == 's':
      return self._Z_source
    elif self._ri == None:
      return None
    else:
      raise Exception(f"Invalid Input Resistance for MOSFET {self._NAME}")
  @ri.setter
  def ri(self,v):
    '''
    Must be 'd', 's', or'g'. 
    Defines what node has the input resistance. 
    '''
    D("", self._DEBUG)
    self._ri = v
    if self._ri != 'd' or self._ri != 's' or self._ri != 'g':
      raise Exception(f"Invalid Input Resistance for MOSFET {self._NAME}")

  @property
  def ro(self):
    '''
    Returns the output resistance along the defined path. 
    '''
    D("", self._DEBUG)
    if self._ro == 'd':
      return self._ZD
    elif self._ro == 'g':
      return self._ZG
    elif self._ro == 's':
      return self._ZS
    elif self._ro == None:
      return None
    else:
      raise Exception(f"Invalid Output Resistance for MOSFET {self._NAME}") 
  @ro.setter
  def ro(self, v):
    '''
    Must be 'd', 's', or'g'. 
    Defines what value is the output resistance. 
    '''
    D("", self._DEBUG)
    self._ro = v
    if self._ro != 'd' and self._ri != 's' and self._ri != 'g':
      raise Exception(f"Invalid Output Resistance for MOSFET {self._NAME}")




# BJT Function
VT=25E-3

class BJT:
  def __init__(self, NAME, TYPE,  I=None, O=None, MIN_VCE=0.2, b=100, VBE=0.7):
    self._NAME = NAME
    self._TYPE = TYPE
    self._GOOD = False
    # Node Voltages
    self._VE = None
    self._VC = None
    self._VB = None
    self._VCE = None
    self._VBE = VBE
    # Node Currents
    self._IE = None
    self._IC = None
    self._IB = None
    # Comp Current
    self._CC = None
    self._MIN_VCE = MIN_VCE
    # Model Values
    self._GM = None
    self._RE = None
    self._RPI = None
    # AC Analysis Values
    self._RI = I
    self._RO = O
    self._Z_emitter = None
    self._Z_base = None
    self._Z_collector = float('inf')
    self._ZE = None
    self._ZB = None
    self._ZC = None
    # other
    self.r = 6
    self._a = b / (b+1)
    self._b = b

  # --- --- ---
  def PrintAttributes(self):
    print(f" --- --- ---\n"
          f"Attributes of Transister {self._NAME}"
          f"\n   Fundamental Attributes: "
          f"\n\tType = {self._TYPE}"
          f"\n\tα = {round(self._a, 4)}"
          f"\n\tβ = {self._b}"
          f"\n   DC Attributes: "
          f"\n\tVE = {PrintUnits(self._VE, 'V')}"
          f" [-] VB = {PrintUnits(self._VB, 'V')}"
          f"\n\tVC = {PrintUnits(self._VC, 'V')}"
          f"\n\tIE = {PrintUnits(self._IE, 'A', print_prefix='m')}"
          f" [-] IC = {PrintUnits(self._IC, 'A', print_prefix='m')}"
          f" [-] IB = {PrintUnits(self._IB, 'A', print_prefix='m')}"
          f"\n   Internal Attributes: "
          f"\n\tgm = {PrintUnits(self._GM, 'Omega', RoundingDigits=2)}"
          f" [-] re = {PrintUnits(self._RE, 'Omega', RoundingDigits=2)}"
          f" [-] rπ = {PrintUnits(self._RPI, 'Omega', RoundingDigits=2)}"
          f"\n   AC Attributes: "
          f"\n\tri = {PrintUnits(self.ri, 'Omega', RoundingDigits=0)}"
          f" [-] ro = {PrintUnits(self.ro, 'Omega', RoundingDigits=0)}"
          f"\n\tZemitter = {PrintUnits(self._Z_emitter, 'Omega', RoundingDigits=0)}"
          f" [-] Zbase = {PrintUnits(self._Z_base, 'Omega', RoundingDigits=0)}"
          f" [-] Zcollector = {PrintUnits(self._Z_collector, 'Omega', RoundingDigits=0)}"
          f"\n\tZE = {PrintUnits(self._ZE, 'Omega', RoundingDigits=0)}"
          f" [-] ZB = {PrintUnits(self.ZB, 'Omega', RoundingDigits=0)}"
          f" [-] ZC = {PrintUnits(self._ZC, 'Omega', RoundingDigits=0)}"
          f"\n   Assumption Values: "
          f"\n\tVCE = {PrintUnits(self._VCE, 'V')}"
          f"\n\tCC = {PrintUnits(self._CC, 'A', print_prefix='m')}"
          f"\n   Verify Assumptions: ")
    VoltageCheck = self.Check_VCE()
    CurrentCheck = self.CheckBaseCurrent()
  # --- --- ---
  def SILENT_CheckBaseCurrent(self):
    if (self._IB != None and self._CC != None):
      Comparitor = abs((self._CC)/(self._IB))
      if (Comparitor >= 10):
        return True
      else:
        return False

  def CheckBaseCurrent(self):
    if (self._IB != None and self._CC != None):
      Comparitor = abs((self._CC)/(self._IB))
      if (Comparitor >= 9.9):
        print(f"Transister {self._NAME} has negligable base Current:\n\t"
              f"{self._NAME}.IB = {PrintUnits(self._IB, 'A', print_prefix='m')}\n\t"
              f"{PrintUnits(self._IB, 'A', print_prefix='m')} << {PrintUnits(self._CC, 'A', print_prefix='m', RoundingDigits=2)}")
        return True
      elif (Comparitor < 9.9):
        print(f"Transister {self._NAME} has non-negligable base Current:\n\t"
              f"{self._NAME}.IB = {PrintUnits(self._IB, 'A', print_prefix='m')} mA\n\t"
              f"{PrintUnits(self._IB, 'A', print_prefix='m')} NOT << {PrintUnits(self._CC, 'A', print_prefix='m', RoundingDigits=2)}")
        return False
      else:
        print("Unknown Error (#1) in CheckBaseCurrent.")
    '''
    elif self._IB == None and self._CC == None:
      print(f"{self._NAME} Comparison Current and Base Current are unknown")
    elif self._IB == None:
      print(f"{self._NAME} Base Current is unknown.")
    elif self._CC == None:
      print(f"{self._NAME} Comparison Current is unknown.")
    else:
      print("Unknown Error (#2) in CheckBaseCurrent.")
    '''
    # print("Unknown Error in CheckBaseCurrent.")

  # ---
  def SILENT_CheckVCE(self):
    if (self._VCE != None and self._VCE > self._MIN_VCE):
      return True
    else:
      return False

  def Check_VCE(self):
    if (self._VCE != None and self._VCE > self._MIN_VCE):
      print(f"Transister {self._NAME} is currently active.\n\t"
            f"{self._NAME}._VCE = {PrintUnits(self._VCE, 'V')}\n\t"
            f"{PrintUnits(self._VCE, 'V')} > {self._MIN_VCE}")
      return True
    if (self._VCE != None and self._VCE < self._MIN_VCE):
      print(f"Transister {self._NAME} is currently inactive.\n\t"
            f"{self._NAME}._VCE = {PrintUnits(self._VCE, 'V')}\n\t"
            f"{PrintUnits(self._VCE, 'V')} < {self._MIN_VCE}")
      return False
  # --- --- ---
  def CheckAssumptions(self):
    SufficientVCE = self.SILENT_CheckVCE()
    NegligableCurrent = self.SILENT_CheckBaseCurrent()
    if (SufficientVCE and NegligableCurrent):
      self._GOOD = True
  # --- --- ---

  @property
  def isGood(self):
    return self._GOOD

  # --- --- ---

  @property
  def Type(self):
    return self._TYPE

  @Type.setter
  def Type(self, Type):
    self._TYPE = Type
    if self._TYPE == "NPN" and (self.VE == None or self.VB == None):
      if self._VE != None:
        self._VB = self._VE + self._VBE
      elif self._VB != None:
        self._VE = self.VB - self._VBE
    elif self._TYPE == "PNP" and (self._VE == None or self._VB == None):
      if self._VE != None:
        self._VB = self.VE - self._VBE
      elif self._VB != None:
        self._VE = self._VB + self._VBE

  # --- --- ---

  @property
  def VE(self):
    return self._VE
  @VE.setter
  def VE(self, v):
    self._VE = v
    if self._TYPE == None:
      print("Transister Type Undefined")
    elif self._TYPE == "NPN" and self._VB == None:
      self._VB = self._VE + self._VBE
    elif self._TYPE == "PNP" and self._VB == None:
      self._VB = self.VE - self._VBE
    if self._VC != None:
      self._VCE = abs(self._VC - self._VE)
    self.CheckAssumptions()
  # --- --- ---

  @property
  def VC(self):
    return self._VC
  @VC.setter
  def VC(self, v):
    self._VC = v
    if self._VE != None:
      self._VCE = abs(self._VC - self._VE)
    self.CheckAssumptions()
  # --- --- ---

  @property
  def VB(self):
    return self._VB
  @VB.setter
  def VB(self, v):
    self._VB = v
    if self._TYPE == None:
      print("Transister Type Undefined")
    elif self._TYPE == "NPN" and self._VE == None:
      self._VE = self._VB - self._VBE
    elif self._TYPE == "PNP" and self._VE == None:
      self._VE = self._VB + self._VBE
    if self._VE != None and self._VC != None:
      self._VCE = abs(self._VC - self._VE)
    self.CheckAssumptions()
  # --- --- ---

  @property
  def IE(self):
    return self._IE
  @IE.setter
  def IE(self, v):
    self._IE = v
    self._IC = v
    self._IB = v / self._b
    self._GM = v / VT
    self._RE = self._a * (VT / v)
    self._RPI = self._b * (VT / v)
    self.CheckAssumptions()
    self.CALC_ZEmitter()
    self.CALC_ZBase()
  # --- --- ---

  @property
  def IC(self):
    return self._IC
  @IC.setter
  def IC(self, v):
    self.IE = v

  # --- --- ---

  @property
  def IB(self):
    return self._IB

  @IB.setter
  def IB(self, v):
    self.IE = (v * self._b)

  # --- --- ---

  @property
  def CC(self):
    return self._CC

  @CC.setter
  def CC(self, v):
    self._CC = v

  # --- --- ---

  @property
  def re(self):
    return self._RE

  @re.setter
  def re(self, v):
    self._RE = v
    self.CALC_ZBase()

  # --- --- ---

  @property
  def rpi(self):
    return self._RPI

  @re.setter
  def rpi(self, v):
    self._RPI = v

  # --- ---

  @property
  def ri(self):
    if self._RI == 'E':
      return self._Z_emitter
    elif self._RI == 'B':
      return self._Z_base
    elif self._RI == 'C':
      return self._Z_collector
    else:
      return None
  @ri.setter
  def ri(self, v):
    if self._RI == None:
      self._RI = v
    else:
      print(f"Value of Input Resistance on transister {self._NAME} is "
            "already set. Reset is Not Possible.")

  # --- ---

  @property
  def ro(self):
    if self._RO == 'E':
      return self._ZE
    elif self._RO == 'B':
      return self._ZB
    elif self._RO == 'C':
      return self._ZC
    else:
      return None

  @ro.setter
  def ro(self, v):
    if self._RO == None:
      self._RO = v
    elif self._RO == 'E':
      self._ZE = v
    elif self._RO == 'B':
      self._ZB = v
    elif self._RO == 'C':
      self._ZC = v
    else:
      print(f"Invalid Input Resistance on transister {self._NAME}")

  # --- ---

  @property
  def Z_emitter(self):
    return self._Z_emitter

  @Z_emitter.setter
  def Z_emitter(self, v=None, auto=True):
    self._Z_emitter = v

  def CALC_ZEmitter(self):
    if self._RE != None and self._ZB != None:
      self._Z_emitter = self._RE + ((self._ZB)/(self.β + 1))

  # --- ---

  @property
  def Z_base(self):
    return self._Z_base

  @Z_base.setter
  def Z_base(self, v):
    self._Z_base = v

  def CALC_ZBase(self):
    if self.ZE != None and self._RE != None:
      self._Z_base = (self.β + 1) * (self.re + self.ZE)
  # --- ---

  @property
  def Z_collector(self):
    return self._Z_collector

  # --- ---

  @property
  def ZE(self):
    return self._ZE

  @ZE.setter
  def ZE(self, v):
    self._ZE = v
    self.CALC_ZBase()

  # --- ---

  @property
  def ZB(self):
    return self._ZB

  @ZB.setter
  def ZB(self, v):
    self._ZB = v
    self.CALC_ZEmitter()

  # --- ---

  @property
  def ZC(self):
    return self._ZC

  @ZC.setter
  def ZC(self, v):
    self._ZC = v

  # --- ---

  @property
  def β(self):
    return self._b

  @β.setter
  def β(self, v):
    self._b = v
    self._a = v / (v + 1)

  # --- ---

  @property
  def _α(self):
    return self._a


# Resistor Class:
#@title Resistor Class
# ----- ----- -----
class R:
  # Initializer
  def __init__(self, NAME, Z=None, I=None, V=None):
    self._NAME = NAME
    self._TYPE = "R"
    # Basic Values
    self._Z = Z
    self._I = I
    self._V = V
    self._ERROR = None
    self.OhmsLaw()

  # --- --- --- Check Resistor Information --- --- ---

  def PrintAttributes(self):
    print(f"Values of Resistor {self._NAME}:"
          f"   R = {PrintUnits(self._Z, 'Omega', 'SCI')}"
          f" [-] I = {PrintUnits(self._I, 'A', 'SCI')}"
          f" [-] VR = {PrintUnits(self._V, 'V', 'SCI')}")

  # --- --- --- Calculate Values --- --- ---
  def OhmsLaw(self):
    if self._Z != None and self._Z != 0 and self._I != None and self._I != 0:
      self._V = self._Z * self._I
    elif self._V != None and self._V != 0 and self._I != None and self._I != 0:
      self._ = self._V / self._I
    elif self._V != None and self._V != 0 and self._Z != None and self._Z != 0:
      self._I = self._V / self._Z

  # --- --- --- Define Value Functions --- --- ---

  @property
  def TYPE(self):
    return self._TYPE


  @property
  def Z(self):
    return self._Z
  @Z.setter
  def Z(self, v):
    if v > 0:
      self._Z = v
      self.OhmsLaw()
    else:
      self._ERROR = True
  
  @property
  def R(self):
    return self._Z
  @R.setter
  def R(self, v):
    if v > 0:
      self._Z = v
      self.OhmsLaw()
    else:
      self._ERROR = True
  

  @property
  def I(self):
    return self._I

  @I.setter
  def I(self, v):
    self._I = v
    self.OhmsLaw()


  @property
  def V(self):
    return self._V

  @V.setter
  def V(self, v):
    self._V = abs(v)
    self.OhmsLaw()


# Capacitor Class:
#@title Capitor Class
# ----- ----- -----
class C:
  # Initializer
  def __init__(self, NAME, Z=None, V=None):
    self._NAME = NAME
    self._TYPE = "C"
    # Basic Values
    self._Z = R
    self._V = V

  # --- --- --- Check Resistor Information --- --- ---

  def PrintAttributes(self, FaradMultiple):
    print(f"Values of Resistor {self._NAME}:"
          f"   C = {PrintUnits(self._Z, 'F',FaradMultiple)}"
          f" [-] VR = {PrintUnits(self._V, 'V')}")

  # --- --- --- Define Value Functions --- --- ---

  @property
  def TYPE(self):
    return self._TYPE

  @property
  def Z(self):
    return self._Z

  @Z.setter
  def Z(self, v):
    if v > 0:
      self._Z = v

  @property
  def V(self):
    return self._V

  @V.setter
  def V(self, v):
    self._V = abs(v)



# Voltage Divider Function:
# @title Voltage Divider


class VD:
  def __init__(self, VH=None, VL=None, ZH=None, ZL=None):
    # Define Divider Values
    self._VH = VH
    self._VL = VL
    self._ZH = ZH
    self._ZL = ZL
    self._TYPE = None
    self._NV = None
    # Define Outcome Values:
    self._GAIN = None
    self._FREQUENCY = None
    self._PASSTYPE = None

    # Calculations:
    self.DividerType()

  # --- --- --- Print Divider Attributes --- --- ---

  def PrintAttributes(self):
    print(f"This is a {self._TYPE} Divider\n"
          f"")
  # --- --- --- Check for Divider Type --- --- ---

  def DividerType(self):
    if (self._ZH.TYPE == 'R' and self._ZL.TYPE == 'R') or (self._ZH.TYPE == 'C' and self._ZL.TYPE == 'C'):
      self._TYPE = 'DC'
      self.CALC_DC()
    elif self._ZH.TYPE == 'R' and self._ZL.TYPE == 'C':
      self._TYPE = 'RC'
      self.CALC_RC()
    elif self._ZH.TYPE == 'C' and self._ZL.TYPE == 'R':
      self._TYPE = 'CR'
      self.CALC_CR()

  # --- --- --- Calculate DC Divided Voltage --- --- ---

  def CALC_DC(self):
    self._ZH.V = ((self._ZH.Z)/(self._ZH.Z + self._ZL.Z)) * abs(self._VH - self._VL)
    self._ZL.V = ((self._ZL.Z)/(self._ZH.Z + self._ZL.Z)) * abs(self._VH - self._VL)
    self._NV = self._VH - self._ZH.V
    self._GAIN = (self._ZH.Z)/(self._ZH.Z + self._ZL.Z)
    self._FREQUENCY = 0
    self._PASSTYPE = 0

  # --- --- ---

  def CALC_RC(self):
    self._GAIN = 1
    self._FREQUENCY = (self._ZH.Z * self._ZL.Z) ** -1
    self._PASSTYPE = -1

  def CALC_CR(self):
    self._GAIN = 1
    self._FREQUENCY = (self._ZH.Z * self._ZL.Z) ** -1
    self._PASSTYPE = 1
