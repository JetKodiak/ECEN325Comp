# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Debug Message
def D(Message, Toggle=False):
  if Toggle:
    print(Message)

# Unit Manipulation
def PrintUnits(v, v_unit, print_prefix=None, v_prefix=None, RoundingDigits=2,debug=False):
  units = ["V", "A", 'Omega', "Ω", "F", '(A/V^2)']
  Prefix = ["E", "P", "T", "G", "M", "k", "h", "da", "",
            "d", "c", "m", "u", "n", "p", "f", "a", 'SCI']
  PrefixValue = [1E18, 1E15, 1E12, 1E9, 1E6, 1E3, 1E2, 1E1, 1,
                  1E-1, 1E-2, 1E-3, 1E-6, 1E-9, 1E-12, 1E-15, 1E-18]

  Standard_Prefixs = ["E", "P", "T", "G", "M", "k", "", "m", "u", "n", "p", "f", "a"]
  Standard_PrefixValue = [1E18, 1E15, 1E12, 1E9, 1E6, 1E3, 1, 0.001, 1E-6, 1E-9, 1E-12, 1E-15, 1E-18]
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

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Operations:
def parallel(ListOfValues):
  InverseSum = 0
  for i in ListOfValues:
    InverseSum += (i ** -1)
  return (InverseSum ** -1)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Resistor Class:
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
    if self._Z is not None and self._Z is not 0 and self._I is not None and self._I is not 0:
      self._V = self._Z * self._I
    elif self._V is not None and self._V is not 0 and self._I is not None and self._I is not 0:
      self._Z = self._V / self._I
    elif self._V is not None and self._V is not 0 and self._Z is not None and self._Z is not 0:
      self._I = self._V / self._Z

  # --- --- --- Define Value Functions --- --- ---

  @property
  def TYPE(self):
    return self._TYPE


  @property
  def Z(self):
    self.OhmsLaw()
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
    self.OhmsLaw()
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
    self.OhmsLaw()
    return self._I

  @I.setter
  def I(self, v):
    self._I = v
    self.OhmsLaw()


  @property
  def V(self):
    self.OhmsLaw()
    return self._V

  @V.setter
  def V(self, v):
    self._V = abs(v)
    self.OhmsLaw()

# Capacitor Class:
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

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BJT Function
class BJT:
  def __init__(self, NAME, TYPE,  I=None, O=None, MIN_VCE=0.2, b=100, VBE=0.7):
    VT=25E-3
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
    self._GM = v / self.VT
    self._RE = self._a * (self.VT / v)
    self._RPI = self._b * (self.VT / v)
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

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MOSFET Simulator  - Vo (Base)
'''
This version ONLY includes variables and a PrintAttributes() function to print all the values in the MOSFET. 
'''
# @title MOS0(NAME, TYPE, Vt, BETA=None, kp=None, W=None, L=None)
class BASEMOS:
  '''
  This is the basic MOS simulator. it ONLY holds the various internal values.
  It does NOT calculate any of the other internal values.
  '''
  def __init__(self, NAME, TYPE, Vt, NHIGH=None, BETA=None, kp=None, W=None, L=None, debug=False):
    self._NAME = NAME
    self._TYPE = TYPE
    self._DEBUG = False
    self._NHIGH = None

    # Define Internal Values:
    self._Vt = Vt
    self._BETA = BETA
    self._kp = kp
    self._W = W
    self._L = L

    # Define the Node Voltages
    self._VS = None
    self._VG = None
    self._VD = None
    # Define the Voltage Differences
    self._VGS = None
    self._VDS = None
    self._VGD = None
    self._Vov = None
    # Define the drain current
    self._ID = None

    # Define the AC Resistances:
    self._gm = None

    self._ZGATE = float('inf')
    self._ZDRAIN = float('inf')
    self._ZSOURCE = None

    self._ZG = None
    self._ZD = None
    self._ZS = None

    # Complex Objects
    self._PVM = None
  # ----- ----- ----- ----- -----


  def PrintAttributes(self, DC=True, AC=True, Region=True):
    a = '   '
    print(f" ----- ----- ----- \n"
          f"Metrics for MOS {self._NAME}\n"
          f"{a}Type = {self._TYPE} [-] Vt = {PrintUnits(self._Vt, 'V', 'SCI')} [-] β = {PrintUnits(self.β, '(A/V^2)', 'SCI')}")
    if DC:
      print(f" --- ----- --- \n"
            f"DC Solutions:\n"
            f"{a}VD = {PrintUnits(self._VD, 'V', 'SCI')} [-] VS = {PrintUnits(self._VS, 'V', 'SCI')} [-] VG = {PrintUnits(self._VG, 'V', 'SCI')}\n"
            f"{a}VDS = {PrintUnits(self.VDS, 'V', 'SCI')} [-] VGS = {PrintUnits(self.VGS, 'V', 'SCI')} [-] VDG = {PrintUnits(self.VDG, 'V', 'SCI')}\n"
            f"{a}Vov = {PrintUnits(self._Vov, 'V', 'SCI')}\n{a}ID = {PrintUnits(self.I, 'A', 'SCI')}")
    if AC:
      print(f" --- ----- --- \n"
            f"AC Solutions:\n"
            f"{a}gm = {PrintUnits(self._gm, 'Omega', 'SCI')}\n"
            f"{a}ZS = {PrintUnits(self._ZS, 'Omega', 'SCI', debug=self._DEBUG)} [-] ZG = {PrintUnits(self._ZG, 'Omega', 'SCI', debug=self._DEBUG)} [-] ZD = {PrintUnits(self._ZD, 'Omega', 'SCI', debug=self._DEBUG)}\n"
            f"{a}ZSource = {PrintUnits(self._ZSOURCE, 'Omega', 'SCI', debug=self._DEBUG)} [-] ZGate = {PrintUnits(self._ZGATE, 'Omega', 'SCI', debug=self._DEBUG)} [-] ZDrain = {PrintUnits(self._ZDRAIN, 'Omega', 'SCI', debug=self._DEBUG)}")
    if Region:
      self.PrintRegion()
    print(f" ----- ----- -----")

  # Check Region:
  def PrintRegion(self):
    a = '   '
    print(f" --- ----- --- ")
    print(f"Checking region of MOSFET {self._NAME}:")
    if self.VGS is None or self.Vt is None:
      # Values are unknown.
      if self.VGS is None and self._Vt is None:
        print(f"{a}Region is unknown: VGS and Vt are unknown.")
      elif self.VGS is None:
        print(f"{a}Region is unknown: VGS is unknown.")
      elif self.Vt is None:
        print(f"{a}Region is unknown: Vt is unknown.")

    elif self.VGS > self.Vt:
      #MOSFET NOT in Cutoff
      if self.VDS is None or self.Vov is None:
        if self.VDS is None and self.Vov is None:
          print(f"{a}Region is unknown: VGS and Vov are Unknown.")
        elif self.VDS is None:
          print(f"{a}Region is unknown: VGS is Unknown.")
        elif self.Vov is None:
          print(f"{a}Region is unknown: Vov is Unknown.")
      elif self.VDS > self.Vov:
        print(f"{a}VDS = {PrintUnits(self.VDS, 'V', 'SCI')}, Vov = {PrintUnits(self.Vov, 'V', 'SCI')}\n"
              f"{a}{PrintUnits(self.VDS, 'V', 'SCI')} > {PrintUnits(self.Vov, 'V', 'SCI')}\n"
              f"{a}MOS is in the active region.")
      elif self.VDS < self.Vov:
        print(f"{a}VDS = {PrintUnits(self.VDS, 'V', 'SCI')}, Vov = {PrintUnits(self.Vov, 'V', 'SCI')}\n"
              f"{a}{PrintUnits(self.VDS, 'V', 'SCI')} < {PrintUnits(self.Vov, 'V', 'SCI')}\n"
              f"{a}MOS is in the Triode Region.")

    elif self.VGS < self.Vt:
      #MOSFET in Cutoff
      print(f"{a}VGS = {PrintUnits(self.VGS, 'V', 'SCI')}, Vt = {PrintUnits(self.Vt, 'V', 'SCI')}")
      print(f"{a}VGS < Vt")
      print(f"{a}MOSFET {self._NAME} is Cutoff. ")



  # ----- ----- ----- ----- -----
  # Define Properties and Setters
  # PROP: VT
  @property
  def Vt(self):
    return self._Vt
  @Vt.setter
  def Vt(self, v):
    self._Vt = v

  # PROP: BETA
  @property
  def BETA(self):
    return self._BETA
  @property
  def β(self):
    return self._BETA
  @property
  def MOSBETA(self):
    return self._BETA
  @BETA.setter
  def BETA(self, v):
    self._BETA = v
  @β.setter
  def β(self, v):
    self._BETA = v
  @MOSBETA.setter
  def MOSBETA(self, v):
    self._BETA = v

  # PROP: VS
  @property
  def VS(self):
    return self._VS
  @VS.setter
  def VS(self, v):
    self._VS = v

  # PROP: VG
  @property
  def VG(self):
    return self._VG
  @VG.setter
  def VG(self, v):
    self._VG = v

  # PROP: VD
  @property
  def VD(self):
    return self._VD
  @VD.setter
  def VD(self, v):
    self._VD = v

  # PROP: VGS
  @property
  def VGS(self):
    return self._VGS
  @property
  def VSG(self):
    return self._VGS
  @VGS.setter
  def VGS(self, v):
    self._VGS = abs(v)
  @VSG.setter
  def VSG(self, v):
    self._VGS = abs(v)

  # PROP: VDS
  @property
  def VDS(self):
    return self._VDS
  @property
  def VSD(self):
    return self._VDS
  @VDS.setter
  def VDS(self, v):
    self._VDS = abs(v)
  @VSD.setter
  def VSD(self, v):
    self._VDS = abs(v)

  # PROP: VGD
  @property
  def VGD(self):
    return self._VGD
  @property
  def VDG(self):
    return self._VGD
  @VGD.setter
  def VGD(self, v):
    self._VGD = abs(v)
  @VDG.setter
  def VDG(self, v):
    self._VGD = abs(v)

  # PROP: Vov
  @property
  def Vov(self):
    return self._Vov
  @property
  def VOV(self):
    return self._Vov
  @property
  def vov(self):
    return self._Vov
  @Vov.setter
  def Vov(self, v):
    self._Vov = v
  @VOV.setter
  def VOV(self, v):
    self._Vov = v
  @vov.setter
  def vov(self, v):
    self._Vov = v

  # PROP: ID
  @property
  def ID(self):
    return self._ID
  @property
  def I(self):
    return self._ID
  @ID.setter
  def ID(self, v):
    self._ID = v
  @I.setter
  def I(self, v):
    self._ID = v

  # PROP: gm
  @property
  def gm(self):
    return self._gm
  @property
  def GM(self):
    return self._GM
  @gm.setter
  def gm(self, v):
    self._gm = v
  @GM.setter
  def GM(self, v):
    self._gm = v

  # PROP: ZGATE
  @property
  def ZGate(self):
    return self._ZGATE
  @property
  def ZGATE(self):
    return self._ZGATE

  # PROP: ZDRAIN
  @property
  def ZDrain(self):
    return self._ZDRAIN
  @property
  def ZDRAIN(self):
    return self._ZDRAIN

  # PROP: ZSOURCE
  @property
  def ZSource(self):
    return self._ZSOURCE
  @property
  def ZSOURCE(self):
    return self._ZSOURCE
  @ZSource.setter
  def ZSource(self, v):
    self._ZSOURCE = v
  @ZSOURCE.setter
  def ZSOURCE(self, v):
    self._ZSOURCE = v

  # PROP: ZG
  @property
  def ZG(self):
    return self._ZG
  @ZG.setter
  def ZG(self, v):
    self._ZG = v

  # PROP: ZD
  @property
  def ZD(self):
    return self._ZD
  @ZG.setter
  def ZD(self, v):
    self._ZD = v

  # PROP: ZS
  @property
  def ZS(self):
    return self._ZS
  @ZS.setter
  def ZS(self, v):
    self._ZS = v

  # ----- ----- ----- ----- -----
# MOSFET Simulator - V1
'''
This version includes a function for calculating VSG when in parallel with another voltage. 
'''
class MOS1(BASEMOS):
  def __init__(self, NAME, TYPE, Vt, BETA=None, kp=None, W=None, L=None, debug=False):
    super().__init__(NAME, TYPE, Vt, BETA, kp, W, L, debug)
    self.CALC_β()
  # ----- ----- ----- ----- -----
    

  # ----- ----- ----- ----- -----
  # Calculate Beta from k, w, and l.
  def CALC_β(self):
    if self.β is None and self._kp is not None and self._W is not None and self._L is not None:
      self.β = self._kp * (self._W / self._L)

  # Find Voltage Drop between nodes.
  def CALC_VΔ(self):
    UPDATE=True
    while UPDATE:
      if self.VDS is None and self.VD is not None and self.VS is not None:
        self.VDS = self.VD - self.VS
        UPDATE = True
      elif self.VGS is None and self.VG is not None and self.VS is not None:
        self.VGS = self.VG - self.VS
        self.Vov = self.VGS - self.Vt
        self.ID = self.β * (1/2) * (self.Vov ** 2)
        UPDATE = True
      elif self.VGD is None and self.VG is not None and self.VD is not None:
        self.VGD = self.VG - self.VD
        UPDATE = True
      else:
        UPDATE = False
  def CALC_NV(self):
    UPDATE=True
    if self._NHigh == 's':
      while (UPDATE):
        if self.VS is None and self.VDS is not None and self.VD is not None:
          UPDATE=True
        elif self.VS is None and self.VGS is not None and self.VG is not None:
          UPDATE=True
        elif self.VG is None and self.VGD is not None and self.VD is not None:
          UPDATE=True
        elif self.VG is None and self.VGS is not None and self.VS is not None:
          UPDATE=True
        elif self.VD is None and self.VDS is not None and self.VS is not None:
          UPDATE=True
        elif self.VD is None and self.VDG is not None and self.VG is not None:
          UPDATE=True
        else:
          UPDATE=False
    elif self._NHIGH == 'd':
      while (UPDATE):
        if self.VS is None and self.VDS is not None and self.VD is not None:
          UPDATE=True
        elif self.VS is None and self.VGS is not None and self.VG is not None:
          UPDATE=True
        elif self.VG is None and self.VGD is not None and self.VD is not None:
          UPDATE=True
        elif self.VG is None and self.VGS is not None and self.VS is not None:
          UPDATE=True
        elif self.VD is None and self.VDS is not None and self.VS is not None:
          UPDATE=True
        elif self.VD is None and self.VDG is not None and self.VG is not None:
          UPDATE=True
        else:
          UPDATE=False
    else:
      return

# MOSFET Simulator - V2
class MOS2(MOS1):
    # ----- ----- ----- ----- -----
  # Get The Quadratic Summation Problems.
  def ParallelVoltageMethod(self, VP, Res, EC=0, PRINT=False):
    '''
    This function takes a parallelvoltage, a resistance (and if it exists, an 
    extra current through the given resistance) and produces the value of
    VGS/VSG on the given MOSFET.
    Inputs: 
    VP (float): This voltage is in parallel with VSG and the voltage across the 
    resistor R.
    Res (float OR R()): This is the resistance that lies before the Source of
    the MOSFET. If it is entered as a Resistance, this function will return a 
    R object which contains the resistance, current, and voltage across the
    resistor. If Res is a R() object then the current and voltage will be
    updated within that object. 
    EC (float): This is any extra current traveling through the resistor R ontop
    of the current from this MOSFET. 
    PRINT (bool): Default = False. If set to True, then the steps of the process
    are printed to console. 
    '''
    if self.VGS != None or self.ID != None:
      raise Exception(f"Parallel Voltage Method is Unnecessary.")
    elif self.Vt == None:
      raise Exception(f"Parallel Voltage Method is either impossible or Unnecessary.")
    ToReturn = None
    if type(Res) == type(float()):
      Res = R("Pre-Resistance", Z=Res)
      ToReturn = Res
    VP = VP - (Res.R * EC)
    # Other Calculations
    A = self.β / 2
    B = ((1 / Res.R) - (self.β * abs(self.Vt)))
    C = (self.β / 2) * (abs(self.Vt) ** 2) - (abs(VP) / Res.R)
    VGSp = ( (-1 * B) + (((B**2) - 4*A*C)) ** (1/2)) / (2 * A)
    VGSm = ( (-1 * B) - (((B**2) - 4*A*C)) ** (1/2)) / (2 * A)
    self.VGS = VGSp
    PVM = (f"ParallelVoltageMethod Used for MOSFET {self._NAME}:\n"
          f"This function uses Quadratic Formula to find the outputs.\n"
          f"Ax^2 + Bx + C = 0 [---] "
          f"A = {A} [-] B = {B} [-] C = {C}\n"
          f"VGSp = {VGSp} [-] VGSm = {VGSm}\n"
          f"VGS = {self.VGS}")
    D(PVM, PRINT)

    # Calculate other related values
    self.Vov = self.VGS - self.Vt
    self.ID = Res.I = (self.β / 2) * (self.Vov ** 2)
    if ToReturn is not None:
      return ToReturn


'''
# --- --- --- 
# MOS1 Test:
kpn = 0.0932174
Vtn = 2.236
M1 = MOS1('M0', 'NMOS', Vt=Vtn, kp=kpn, W=1, L=1)
M1.VG = 5
M1.VS = 1
M1.CALC_VΔ()
M1.PrintAttributes()
'''