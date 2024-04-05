# @title Print w/ Units
units = ["V", "A", 'Omega', "Ω", "F"]
Prefix = ["E", "P", "T", "G", "M", "k", "h", "da", "-",
          "d", "c", "m", "u", "n", "p", "f", "a"]
PrefixValue = [1E18, 1E15, 1E12, 1E9, 1E6, 1E3, 1E2, 1E1, 1,
                1E-1, 1E-2, 1E-3, 1E-6, 1E-9, 1E-12, 1E-15, 1E-18]
def PrintUnits(v, v_unit, print_prefix=None, v_prefix=None, RoundingDigits=4):
  # Major Errors:
  if v == None or v_unit == None:
    return "[Not Defined]"
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
    v_prefix = "-"
  if unit in units and v_prefix in Prefix and print_prefix in Prefix:
    v_index = Prefix.index(v_prefix)
    print_index = Prefix.index(print_prefix)
    vi = v * PrefixValue[v_index]
    vp = vi / PrefixValue[print_index]
    return f"{round(vp, RoundingDigits)} {print_prefix}{unit}"
  # --- --- ---
  # ERror Detection
  if not unit in units:
    print("Invalid Unit")
  elif not v_prefix in Prefix:
    print(f"Invalid Value Prefix: {v_prefix}")
  elif not print_prefix in Prefix:
    print("Invalid print prefix")

# --- --- ---

# @title Operations
def parallel(ListOfValues):
  InverseSum = 0
  for i in ListOfValues:
    InverseSum += (i ** -1)
  return (InverseSum ** -1)

'''
The following function represents a voltage divider where the higher voltage 
is attached to R1, and the lower voltage is attached to R2. V1 is the higher 
voltage and R2 is the lower voltage
'''

# --- --- ---

# @title MOS Class
class MOS:
  def __init__(self, NAME, TYPE, Width=None, Length=None, KPRIME = None, MOSBETA = None, Vtn = None):
    # Define Fundamental Elements
    self._NAME = NAME
    self._TYPE = TYPE
    # Component Values:
    self._Vtn = Vtn
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
    self._VD = None
    self._VG = None
    self._VS = None
    # Define Voltage Differences:
    self._VDS = None
    self._VGS = None
    self._VDG = None
    self._Vov = None
    # Define the current:
    self._ID = None
  
  # --- --- --- Print MOS Metrics --- --- ---

  def PrintAttributes(self):
    a = '\t'
    print(f" --- --- --- \n"
          f"Metrics for MOS {self._NAME}\n"
          f"{a}Type = {self._TYPE} [-] Vtn = {self._Vtn} [-] MOSβ = {self._MOSBETA}\n"
          f"{a}VD = {self._VD} [-] VS = {self._VS} [-] VG = {self._VG}\n"
          f"{a}VDS = {self._VDS} [-] VGS = {self._VGS} [-] VDG = {self._VDG}\n"
          f"{a}Vov = {self._Vov} [-] ID = {self._ID}")
  
  def Active(self):
    self.CheckRegion()
    if self._CUTOFF:
      print(f"VGS = {self._VGS}, Vtn = {self._Vtn}\n"
            f"VGS < Vtn\n"
            f"MOS is in Cutoff Region.")
    elif self._SAT:
      print(f"VDS = {self._VDS}, Vov = {self._Vov}\n"
            f"{self._VDS} > {self._Vov}\n"
            f"MOS is in active region.")
    elif not self._SAT:
      print(f"VDS = {self._VDS}, Vov = {self._Vov}\n"
            f"{self._VDS} < {self._Vov}\n"
            f"MOS is in Triode Region.")
    else: 
      print(f"VDS = {self._VDS}, Vov = {self._Vov}\n"
            f"MOS region is unknown.")

  # --- --- --- Check whether is Active region --- --- ---

  def CheckRegion(self):
    if self._VSG < self._Vtn:
      self._CUTOFF = True
      self._SAT = False
    elif self._VDS != None and self._Vov != None:
      self._SAT = self._VDS > self._Vov
      self._CUTOFF = False
  
  # --- --- --- Calculate Values --- --- ---

  def CALC_MOSBETA(self):
    '_MOSBETA depends on L, W, and kprime'
    if self._MOSBETA == None and self._KPRIME != None and self._W != None and self._L != None:
      self._MOSBETA = self._KPRIME * (self._W / self._L)

  def CALC_Vov(self):
    'Vov depends on VGS and Vtp'
    if self._Vov == None and self._VGS != None and self._Vtp != None:
      self._Vov = abs(self._VGS) - abs(self._Vtn)
    if self._Vov != None and self._ID != None and self._MOSBETA != None:
      self._Vov = ((self._ID * 2)/(self._MOSBETA)) ** (1/2)
    if self._Vov != None:
      self.CheckRegion()

  def CALC_VDS(self):
    'VDS depends on VD and VS'
    if self._VD != None and self._VS != None:
      self._VDS = abs(self._VD - self._VS)
      self.CheckRegion()

  def CALC_VGS(self):
    'VGS depends on VG and VS'
    if self._VG != None and self._VS != None:
      self._VGS = abs(self._VG - self._VS)

  def CALC_VDG(self):
    'VDG depends on VD and VG'
    if self._VD != None and self._VG != None:
      self._VDG = abs(self._VD - self._VG)

  def CALC_ID(self):
    'ID depends on MOSBETA and Vov'
    if self._SAT == None:
      print(f"Saturation is Unknown")
    elif not self._SAT:
      print(f"MOSFET is NOT active. ")
    if self._Vov != None and self._MOSBETA != None and self._SAT:
      self._ID = (self._MOSBETA * (self._Vov ** 2))/(2)

  # --- --- --- Define all the properties and setters --- --- ---
  # --- Define Voltage Functions --- 
  
  @property
  def VD(self):
    return self._VD

  @VD.setter
  def VD(self, v):
    self._VD = v
    # Define Related Values
    self.CALC_VDS()
    self.CALC_VDG()
  
  @property
  def VG(self):
    return self._VG
  @VD.setter
  def VG(self, v):
    self._VG = v
    # Calculate any related values
    self.CALC_VGS()
    self.CALC_VDG()
  
  @property
  def VS(self):
    return self._VS
  @VS.setter
  def VS(self, v):
    self._VS = v
    self.CALC_VDS()
    self.CALC_VGS()
  
  # --- Define Voltage Difference Functions ---

  @property
  def VDS(self):
    return self._VDS
  
  @property
  def VGS(self):
    return self._VGS

  @property
  def VDG(self):
    return self._VDG
  
  # --- Define Current Functions ---

  @property
  def ID(self):
    return self._ID
  @ID.setter
  def ID(self, v):
    self._ID = v
    self.CALC_Vov()

# --- --- ---

#@title BJT Class
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

# --- --- ---

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
          f"   R = {PrintUnits(self._Z, 'Omega')}"
          f" [-] I = {PrintUnits(self._I, 'A', 'm')}"
          f" [-] VR = {PrintUnits(self._V, 'V')}")
  
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

# --- --- ---

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

# --- --- ---

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
    self._NV = self._VL + self._ZH.V
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

# --- --- --
