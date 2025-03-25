#region imports
import sys
from ThermoStateCalc import Ui__frm_StateCalculator
from pyXSteam.XSteam import XSteam
from PyQt5.QtWidgets import QWidget, QApplication
from UnitConversion import UC
from scipy.optimize import fsolve
#endregion

#region class definitions
class thermoSatProps:
    def __init__(self, p=None, t=None, SI=True):
        '''
        This is a class to compute the thermodynamic properties of water as:
        1. saturated liquid
        2. two-phase
        3. saturated vapor
        4. superheated vapor
        :param p: pressure in appropriate units
        :param t: temperature in appropriate units
        :param SI: boolean True=SI units, False = English units
        '''
        #creates a steamTable object from XSteam
        self.steamTable=XSteam(XSteam.UNIT_SYSTEM_MKS)
        if(p is not None):
            self.getSatProps(p, SI)
        elif(t is not None):
            self.getSatProps(self.steamTable.psat_t(t),SI)

    def getSatProps(self, p, SI=True):
        '''
        Retrieve the saturated liquid and saturated vapor properties at a specified pressure.
        :param p:
        :param SI:
        :return:
        '''
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        # compute saturated properties at isobar p
        self.pSat = p
        self.tSat = self.steamTable.tsat_p(p)
        self.vf = self.steamTable.vL_p(p)
        self.vg = self.steamTable.vV_p(p)
        self.hf = self.steamTable.hL_p(p)
        self.hg = self.steamTable.hV_p(p)
        self.uf = self.steamTable.uL_p(p)
        self.ug = self.steamTable.uV_p(p)
        self.sf = self.steamTable.sL_p(p)
        self.sg = self.steamTable.sV_p(p)
        self.vgf = self.vg - self.vf
        self.hgf = self.hg - self.hf
        self.sgf = self.sg - self.sf
        self.ugf = self.ug - self.uf
class thermoState:
    def __init__(self, p=None,t=None,v=None,u=None,h=None,s=None,x=None):
        """
        This is a class I use for storing a thermodynamic state.  Calling setState requires you to specify two
        independent thermodynamic properties.  One ambiguity exists if you specify both psat and tsat.  In that
        case I assume two-phase with x=0.5
        :param p:
        :param t:
        :param v:
        :param u:
        :param h:
        :param s:
        :param x:
        """
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.region = "saturated"
        self.p=p
        self.t=t
        self.v=v
        self.u=u
        self.h=h
        self.s=s
        self.x=x

    def computeProperties(self):
        """
        Assumes p and t are already calculated
        :return:
        """
        if self.region ==  "two-phase":
            #since two-phase, use quality to interpolate overall properties
            self.u = self.steamTable.uL_p(self.p)+self.x*(self.steamTable.uV_p(self.p)-self.steamTable.uL_p(self.p))
            self.h = self.steamTable.hL_p(self.p) + self.x * (self.steamTable.hV_p(self.p) - self.steamTable.hL_p(self.p))
            self.s = self.steamTable.sL_p(self.p) + self.x * (self.steamTable.sV_p(self.p) - self.steamTable.sL_p(self.p))
            self.v = self.steamTable.vL_p(self.p) + self.x * (self.steamTable.vV_p(self.p) - self.steamTable.vL_p(self.p))
        else:
            self.u = self.steamTable.u_pt(self.p, self.t)
            self.h = self.steamTable.h_pt(self.p, self.t)
            self.s = self.steamTable.s_pt(self.p, self.t)
            self.v = self.steamTable.v_pt(self.p, self.t)
            self.x = 1.0 if self.region == "super-heated vapor" else 0.0

    def setState(self, stProp1, stProp2, stPropVal1, stPropVal2, SI=True):
        """
        Calculates the thermodynamic state variables based on specified values.
        I have thermodynamic variables:  P, T, v, h, u, s and x (7 things) from which I am choosing two.
        Possible number of permutations:  7!/5! =42.
        But, order of the two things does not matter, so 42/2=21
        PT, Pv, Ph, Pu, Ps, Px (6)
        Tv, Th, Tu, Ts, Tx (5)
        vh, vu, vs, vx (4)
        hu, hs, hx (3)
        us, ux (2)
        sx (1)
        Total of 21 cases to deal with.  I will attack them in the order shown above
        :return: nothing
        """
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        # Step 1: read which properties are being specified from the combo boxes
        SP = [stProp1, stProp2]
        SP[0] = SP[0].lower()
        SP[1] = SP[1].lower()
        f1 = float(stPropVal1)
        f2 = float(stPropVal2)
        # Step 2: select the proper case from the 21.  Note that PT is the same as TP etc.
        if SP[0] == 'p' or SP[1] == 'p':
            oFlipped = SP[0] != 'p'
            SP1 = SP[0] if oFlipped else SP[1]
            self.p = f1 if not oFlipped else f2
            tSat = self.steamTable.tsat_p(self.p)
            # case 1:  pt or tp
            if SP1 == 't':
                self.t = f2 if not oFlipped else f1
                tSat = round(tSat)  # I will compare at 3 three decimal places
                # compare T to TSat
                if self.t < tSat or self.t > tSat:
                    self.region = "sub-cooled liquid" if self.t < tSat else "super-heated vapor"
                else:  # this is ambiguous since at saturated temperature
                    self.region = "two-phase"
                    self.x = 0.5
            # case 2: pv or vp
            elif SP1 == 'v':
                self.v = f2 if not oFlipped else f1
                vf = round(self.steamTable.vL_p(self.p), 5)
                vg = round(self.steamTable.vV_p(self.p), 3)
                # compare v to vf and vg
                if self.v < vf or self.v > vg:
                    self.region = "sub-cooled liquid" if self.v < vf else "super-heated vapor"
                    # since I can't find properties using v, I will use fsolve to find T
                    dt = 1.0 if self.v > vg else -1.0
                    fn1 = lambda T: self.v - self.steamTable.v_pt(self.p, T)
                    self.t = fsolve(fn1, [tSat + dt])[0]
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.v - self.vf) / (self.vgf)
                    self.t = tSat
            # case 3 pu or up
            elif SP1 == 'u':
                self.u = f2 if not oFlipped else f1
                uf = round(self.steamTable.uL_p(self.p), 5)
                ug = round(self.steamTable.uV_p(self.p), 3)
                ugf = ug-uf
                # compare u to uf and ug
                if self.u < uf or self.u > ug:
                    self.region = "sub-cooled liquid" if self.u < uf else "super-heated vapor"
                    # since I can't find properties using u, I will use fsolve to find T
                    dt = 1.0 if self.u > ug else -1.0
                    fn3 = lambda T: self.u - self.steamTable.u_pt(self.p, T)
                    self.t = fsolve(fn3, [tSat + dt])[0]
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - uf) / (ugf)
                    self.t = tSat
            # case 4 ph or hp
            elif SP1 == 'h':
                self.h = f2 if not oFlipped else f1;
                hf = self.steamTable.hL_p(self.p)
                hg = self.steamTable.hL_p(self.p)
                hgf=hg-hf
                # compare h to hf and hg
                if self.h < hf or self.h > hg:
                    self.region = "sub-cooled liquid" if self.h < hf else "super-heated vapor"
                    self.t = self.steamTable.t_ph(self.p, self.h)
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.h - hf) / (hgf)
                    self.t = tSat
            # case 5 ps or sp
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                sf=self.steamTable.sL_p(self.p)
                sg=self.steamTable.sV_p(self.p)
                sgf=sg-sf
                # compare s to sf and sg
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                    self.t = self.steamTable.t_ps(self.p, self.s)
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.s - sf) / (sgf)
            # case 6 px or xp
            elif SP1 == 'x':
                self.region = "two-phase"
                self.x = f2 if not oFlipped else f1
                self.t = tSat
        elif SP[0] == 't' or SP[1] == 't': # t and another property specified
            oFlipped = SP[0] != 't'
            SP1 = SP[0] if oFlipped else SP[1]
            self.t=f1 if not oFlipped else f2
            pSat = self.steamTable.psat_t(self.t)
            # case 7:  tv or vt
            if SP1 == 'v':
                self.v = f2 if not oFlipped else f1
                vf=self.steamTable.vL_p(pSat)
                vg=self.steamTable.vV_p(pSat)
                vgf=vg-vf
                # compare v to vf and vg
                if self.v < vf or self.v > vg:
                    self.region = "sub-cooled liquid" if self.v < vf else "super-heated vapor"
                    # since I can't find properties using v, I will use fsolve to find P
                    dp = -0.1 if self.v > vg else 0.1
                    fn3 = lambda P: [self.v - self.steamTable.v_pt(P, self.t)]
                    self.p = fsolve(fn3, [pSat + dp])[0]
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.v - vf) / (vgf)
                    self.p = self.pSat
            # case 8:  tu or ut
            elif SP1 == 'u':
                self.u = f2 if not oFlipped else f1
                uf=self.steamTable.uL_p(pSat)
                ug=self.steamTable.uV_p(pSat)
                ugf=ug-uf
                # compare u to uf and ug
                if self.u < uf or self.u > ug:
                    self.region = "sub-cooled liquid" if self.u < uf else "super-heated vapor"
                    # since I can't find properties using u, I will use fsolve to find P
                    dp = 0.1 if self.u > ug else -0.1
                    fn8 = lambda P: self.u - self.steamTable.u_pt(self.t, P)
                    self.p = fsolve(fn8, [pSat + dp])[0]
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - uf) / (ugf)
                    self.p = self.pSat
            # case 9:  th or ht
            elif SP1 == 'h':
                self.h = f2 if not oFlipped else f1
                hf=self.steamTable.hL_p(pSat)
                hg=self.steamTable.hV_p(pSat)
                hgf=hg-hf
                # compare h to hf and hg
                if self.h < hf or self.h > hg:
                    self.region = "sub-cooled liquid" if self.h < hf else "super-heated vapor"
                    self.p = self.steamTable.p_th(self.t, self.h)
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = pSat
                    self.x = (self.h - hf) / (hgf)
            # case 10:  ts or st
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                sf=self.steamTable.sL_p(pSat)
                sg=self.steamTable.sV_p(pSat)
                sgf=sg-sf
                # compare s to sf and sg
                if self.s < self.sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                    self.p = self.steamTable.p_ts(self.t, self.s)
                    # now use P and T
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = pSat
                    self.x = (self.s - sf) / (sgf)
            # case 11:  tx or xt
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                self.p = pSat
        elif SP[0] == 'v' or SP[1] == 'v':
            oFlipped = SP[0] != 'v'
            SP1 = SP[0] if oFlipped else SP[1]
            self.v = f1 if not oFlipped else f2
            # case 12:  vh or hv
            if SP1 == 'h':
                self.h = f2 if not oFlipped else f1
                #find p where h and v match
                def fn12(P):
                    # could be single phase or two-phase, but both v&h have to match at same x
                    hf=self.steamTable.hL_p(P)
                    hg=self.steamTable.hV_p(P)
                    hgf=hg-hf
                    vf=self.steamTable.vL_p(P)
                    vg=self.steamTable.vV_p(P)
                    vgf=vg-vf
                    if self.between(self.h, hf, hg):
                        self.x = (self.h - hf) / hgf
                        return self.v - (vf + self.x * vgf)
                    # could be single phase
                    return self.v - self.steamTable.v_ph(P, self.h)

                self.p = fsolve(fn12, [1.0])[0]
                vf=self.steamTable.vL_p(self.p)
                vg=self.steamTable.vV_p(self.p)
                tsat=self.steamTable.tsat_p(self.p)
                # compare v to vf and vg
                if self.v < vf or self.v > vg:
                    self.region = "sub-cooled liquid" if self.v < vf else "super-heated vapor"
                    dt = -1 if self.v<vf else 1
                    findtgivenv=lambda t: self.v-self.steamTable.v_pt(self.p, t)
                    self.t=fsolve(findtgivenv,[tsat+dt])[0]
                else:  # two-phase
                    self.region = "two-phase"
                    self.t=tsat
                    # first calculate quality
                    self.x = (self.v - self.vf) / (self.vgf)
            # case 13:  vu or uv
            elif SP1 == 'u':
                self.u = f2 if not oFlipped else f1

                # use fsolve to find P&T at this v & u
                def fn13(PT):
                    p,t=PT
                    uf=self.steamTable.uL_p(p)
                    ug=self.steamTable.uV_p(p)
                    ugf=ug-uf
                    vf=self.steamTable.vL_p(p)
                    vg=self.steamTable.vV_p(p)
                    vgf=vg-vf
                    if self.between(self.u, uf,ug):
                        self.t = self.tSat
                        self.x = (self.u - uf) / ugf
                        return [self.v - (vf + self.x * vgf), 0]
                    return [self.v - self.steamTable.v_pt(p, t), self.u - self.steamTable.u_pt(p, t)]

                props = fsolve(fn13, [1, 100])
                self.p = props[0]
                self.t = props[1]
                uf=self.steamTable.uL_p(self.p)
                ug=self.steamTable.uV_p(self.p)
                ugf=ug-uf
                if self.u<uf or self.u>ug:
                    self.region = "sub-cooled" if self.u<uf else "super-heated"
                else:
                    self.region = "two-phase"
                    self.x = (self.u-uf)/ugf
            # case 14:  vs os sv
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                def fn14(PT):
                    p, t = PT
                    sf = self.steamTable.uL_p(p)
                    sg = self.steamTable.uV_p(p)
                    sgf = sg - sf
                    vf = self.steamTable.vL_p(p)
                    vg = self.steamTable.vV_p(p)
                    vgf = vg - vf
                    if self.between(self.s, sf, sg):
                        self.x = (self.s - sf) / sgf
                        return [self.v - vf - self.x * vgf, 0.0]
                    return [self.v - self.steamTable.v_pt(p,t),
                            self.s - self.steamTable.s_pt(p,t)]
                props = fsolve(fn14, [1, self.steamTable.sV_p(1)])
                self.p = props[0]
                self.t = props[1]
                sg=self.steamTable.sV_p(self.p)
                sf=self.steamTable.sL_p(self.p)
                sgf=sg-sf
                # compare s to sf and sg
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - sf) / (sgf)
            # case 15:  vx or xv
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0.0, 1.0)
                self.region = "two-phase"

                def fn15(p):
                    vf=self.steamTable.vL_p(p)
                    vg=self.steamTable.vV_p(p)
                    vgf=vg-vf
                    return self.v - (vf + self.x * vgf)

                props = fsolve(fn15, [1])
                self.p = props[0]
                self.t = self.steamTable.tsat_p(self.p)
        elif SP[0] == 'h' or SP[1] == 'h':
            oFlipped = SP[0] != 'h'
            SP1 = SP[0] if oFlipped else SP[1]
            self.h=f1 if not oFlipped else f2
            # case 16:  hu or uh
            if SP1 == 'u':
                self.u = f2 if not oFlipped else f1
                # use fsolve to find P&T at this h & u
                def fn16(PT):
                    p,t=PT
                    uf=self.steamTable.uL_p(p)
                    ug=self.steamTable.uV_p(p)
                    ugf=ug-uf
                    hf=self.steamTable.hL_p(p)
                    hg=self.steamTable.hV_p(p)
                    hgf=hg-hf
                    if self.between(self.u, uf, ug):
                        self.x = (self.u - self.uf) / self.ugf
                        return [self.h - hf - self.x * hg, 0.0]
                    return [self.h - self.steamTable.h_pt(p,t), self.u - self.steamTable.u_pt(p,t)]

                props = fsolve(fn16, [1, 100])
                self.p = props[0]
                self.t = props[1]
                uf = self.steamTable.uL_p(self.p)
                ug = self.steamTable.uV_p(self.p)
                ugf = ug - uf
                # compare u to uf and ug
                if self.u < uf or self.u > ug:
                    self.region = "sub-cooled liquid" if self.u < uf else "super-heated vapor"
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - uf) / (ugf)
            # case 17:  hs or sh
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                def fn17(PT):
                    p,t=PT
                    sf=self.steamTable.sL_p(p)
                    sg=self.steamTable.sV_p(p)
                    sgf=sg-sf
                    hf=self.steamTable.hL_p(p)
                    hg=self.steamTable.hV_p(p)
                    hgf=hg-hf
                    if self.between(self.s, sf, sg):
                        self.x = (self.s - self.sf) / self.sgf
                        return [self.h - hf - self.x * hg, 0.0]
                    return [self.h - self.steamTable.h_pt(p,t), self.s - self.steamTable.s_pt(p,t)]

                props = fsolve(fn17, [1, 100])
                self.p = props[0]
                self.t = props[1]
                sf = self.steamTable.sL_p(self.p)
                sg = self.steamTable.sV_p(self.p)
                sgf = sg - sf
                # compare s to sf and sg
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.s - sf) / (sgf)
            # case 18:  hx or xh
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0.0, 1.0)
                self.region = "two-phase"

                def fn18(p):
                    hf = self.steamTable.hL_p(self.p)
                    hg = self.steamTable.hV_p(self.p)
                    hgf = hg - hf
                    return self.h - (hf + self.x * hgf)

                props = fsolve(fn18, [1])
                self.p = props[0]
                self.t = self.steamTable.tsat_p(self.p)
        elif SP[0] == 'u' or SP[1] == 'u':
            oFlipped = SP[0] != 'u'
            SP1 = SP[0] if oFlipped else SP[1]
            self.u = f1 if not oFlipped else f2

            # case 19:  us or su
            if SP1 == 's':
                self.s = f2 if not oFlipped else f1
                def fn19(PT):
                    p,t=PT
                    sf=self.steamTable.sL_p(p)
                    sg=self.steamTable.sV_p(p)
                    sgf=sg-sf
                    uf=self.steamTable.uL_p(p)
                    ug=self.steamTable.uV_p(p)
                    ugf=ug-uf
                    if self.between(self.s, sf, sg):
                        self.x = (self.s - self.sf) / self.sgf
                        return [0.0, self.s - sf - self.x * sg]
                    return [self.u - self.steamTable.u_pt(p,t), self.s - self.steamTable.s_pt(p,t)]

                props = fsolve(fn19, [1, 100])
                self.p = props[0]
                self.t = props[1]
                sf = self.steamTable.sL_p(self.p)
                sg = self.steamTable.sV_p(self.p)
                sgf = sg - sf
                # compare s to sf and sg
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.s - sf) / (sgf)
            # case 20:  ux or xu
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x=self.clamp(self.x,0,1)
                self.region = "two-phase"
                def fn20(p):
                    hf=self.steamTable.hL_p(p)
                    hg=self.steamTable.hV_p(p)
                    hgf=hg-hf
                    return self.h - (hf + self.x * hgf)

                props = fsolve(fn20, [1])
                self.p = props[0]
                self.t = self.steamTable.tsat_p(self.p)
        elif SP[0] == 's' or SP[1] == 's':
            oFlipped = SP[0] != 's'
            SP1 = SP[0] if oFlipped else SP[1]
            self.s = f1 if not oFlipped else f2
            # case 21:  sx or xs
            if SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x=self.clamp(self.x,0,1)
                self.region = "two-phase"

                def fn21(p):
                    sf=self.steamTable.sL_p(p)
                    sg=self.steamTable.sV_p(p)
                    sgf = sg - sf
                    return self.s - (sf + self.x * sgf)

                props = fsolve(fn21, [1])
                self.p = props[0]
                self.t = self.steamTable.tsat_p(self.p)
        self.computeProperties()

    def __sub__(self, other):
        delta = thermoState()
        delta.p=self.p-other.p
        delta.t=self.t-other.timeData
        delta.h=self.h-other.h
        delta.u=self.u-other.u
        delta.s=self.s-other.s
        delta.v=self.v-other.v
        return delta
class main_window(QWidget,Ui__frm_StateCalculator):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.SetupSlotsAndSignals()
        self.steamTable=XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.currentUnits='SI'
        self.setUnits()
        self.show()
    def SetupSlotsAndSignals(self):
        """
        I've modified the original to include a state2.  Here I assign slots to GUI actions.
        :return:
        """
        self._rdo_English.clicked.connect(self.setUnits)
        self._rdo_SI.clicked.connect(self.setUnits)
        self._cmb_Property1.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property2.currentIndexChanged.connect(self.setUnits)
        self._pb_Calculate.clicked.connect(self.calculateProperties)

    def setUnits(self):
        """
        This sets the units for the selected specified properties.
        Units for the thermodynamic properties are set upon pushing calculate button.
        :return:
        """
        #set the units system based on selected radio button
        #also, determine if a units change is required
        SI=self._rdo_SI.isChecked()
        newUnits='SI' if SI else 'EN'
        UnitChange = self.currentUnits != newUnits  # compare new units to current units
        self.currentUnits = newUnits

        if SI:
            self.steamTable=XSteam(XSteam.UNIT_SYSTEM_MKS)
            self.l_Units = "m"
            self.p_Units = "bar"
            self.t_Units = "C"
            self.m_Units = "kg"
            self.time_Units = "s"
            self.energy_Units = "W"
            self.u_Units = "kJ/kg"
            self.h_Units = "kJ/kg"
            self.s_Units = "kJ/kg*C"
            self.v_Units = "m^3/kg"
        else:
            self.steamTable=XSteam(XSteam.UNIT_SYSTEM_FLS)
            self.l_Units = "ft"
            self.p_Units = "psi"
            self.t_Units = "F"
            self.m_Units = "lb"
            self.time_Units = "s"
            self.energy_Units = "btu"
            self.u_Units = "btu/lb"
            self.h_Units = "btu/lb"
            self.s_Units = "btu/lb*F"
            self.v_Units = "ft^3/lb"

        #read selected Specified Properties from combo boxes
        SpecifiedProperty1 = self._cmb_Property1.currentText()
        SpecifiedProperty2 = self._cmb_Property2.currentText()
        #read numerical values for selected properties
        SP=[float(self._le_Property1.text()), float(self._le_Property2.text())]

        #region set units labels and convert values if needed for state 1
        #region property 1
        if 'Pressure' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.p_Units)
            if UnitChange:  # note that I only should convert if needed.  Not if I double click on SI or English
                SP[0]=SP[0]*UC.psi_to_bar if SI else SP[0]*UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.t_Units)
            if UnitChange:
                SP[0] = UC.F_to_C(SP[0]) if SI else UC.C_to_F(SP[0])
        elif 'Energy' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.u_Units)
            if UnitChange:
                SP[0]=SP[0]*UC.btuperlb_to_kJperkg if SI else SP[0]*UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.h_Units)
            if UnitChange:
                SP[0] = SP[0] * UC.btuperlb_to_kJperkg if SI else SP[0] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.s_Units)
            if UnitChange:
                SP[0] = SP[0] * UC.btuperlbF_to_kJperkgC if SI else SP[0] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.v_Units)
            if UnitChange:
                SP[0]=SP[0]*UC.ft3perlb_to_m3perkg if SI else SP[0]*UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText("")
        #endregion
        #region property 2
        if 'Pressure' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.p_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.psi_to_bar if SI else SP[1] * UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.t_Units)
            if UnitChange:
                SP[1] = UC.F_to_C(SP[1]) if SI else UC.C_to_F(SP[1])
        elif 'Energy' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.u_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.btuperlb_to_kJperkg if SI else SP[1] * UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.h_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.btuperlb_to_kJperkg if SI else SP[1] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.s_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.btuperlbF_to_kJperkgC if SI else SP[1] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.v_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.ft3perlb_to_m3perkg if SI else SP[1] * UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText("")
        #endregion

        self._le_Property1.setText("{:0.3f}".format(SP[0]))
        self._le_Property2.setText("{:0.3f}".format(SP[1]))
        #endregion

    def clamp(self, x, low, high):
        """
        This clamps a float x between a high and low limit inclusive
        :param x:
        :param low:
        :param high:
        :return:
        """
        if x<low:
            return low
        if x>high:
            return high
        return x

    def between(self, x, low, high):
        """
        Tells if x is between low and high inclusive
        :param x:
        :param low:
        :param high:
        :return:
        """
        if x>=low and x<=high:
            return True
        return False

    def getSatProps_p(self, p):
        """
        Given a pressure, calculate the saturated properties for that isobar
        :param p:
        :return: a thermoSatProps object containing all the saturated properties for that isobar
        """
        return thermoSatProps(p=p)

    def getSatProps_t(self, t):
        """
        Given a temperature, calculate the saturation pressure and then
        calculate all other saturated properties
        :param t:
        :return: a thermoSatProps object containing all the saturation properties at this temperature
        """
        return thermoSatProps(t=t)

    def makeLabel(self, state):
        """
        Given that I have T&P, find the other properties and make the label
        :return:
        """
        stProps = "Region = {:}".format(state.region)
        stProps += "\nPressure = {:0.3f} ({:})".format(state.p, self.p_Units)
        stProps += "\nTemperature = {:0.3f} ({:}) ".format(state.t, self.t_Units)
        stProps += "\nInternal Energy = {:0.3f} ({:})".format(state.u, self.u_Units)
        stProps += "\nEnthalpy = {:0.3f} ({:})".format(state.h, self.h_Units)
        stProps += "\nEntropy = {:0.3f} ({:})".format(state.s, self.s_Units)
        stProps += "\nSpecific Volume = {:0.3f} ({:})".format(state.v, self.v_Units)
        stProps += "\nQuality = {:0.3f}".format(state.x)
        return stProps
    def makeDeltaLabel(self, state1,state2):
        """
        calculatet the change in thermodynamic properties and put them in the form of a label
        :param state1:
        :param state2:
        :return:
        """
        stDelta="Property change:"
        stDelta+="\nT2-T1 = {:0.3f} {:}".format(state2.timeData - state1.timeData, self.t_Units)
        stDelta+="\nP2-P1 = {:0.3f} {:}".format(state2.p-state1.p, self.p_Units)
        stDelta+="\nh2-h1 = {:0.3f} {:}".format(state2.h-state1.h, self.h_Units)
        stDelta+="\nu2-u1 = {:0.3f} {:}".format(state2.u-state1.u, self.u_Units)
        stDelta+="\ns2-s1 = {:0.3f} {:}".format(state2.s-state1.s, self.s_Units)
        stDelta+="\nv2-v1 = {:0.3f} {:}".format(state2.v-state1.v, self.v_Units)
        return stDelta

    def calculateProperties(self):
        """
        Calculates the thermodynamic state variables based on specified values.
        I have thermodynamic variables:  P, T, v, h, u, s and x (7 things) from which I am choosing two.
        Possible number of permutations:  7!/5! =42.
        But, order of the two things does not matter, so 42/2=21
        PT, Pv, Ph, Pu, Ps, Px (6)
        Tv, Th, Tu, Ts, Tx (5)
        vh, vu, vs, vx (4)
        hu, hs, hx (3)
        us, ux (2)
        sx (1)
        Total of 21 cases to deal with.  I will attack them in the order shown above
        :return: nothing
        """
        self.state1 = thermoState()
        self.state2 = thermoState()

        SP=[self._cmb_Property1.currentText()[-2:-1].lower(), self._cmb_Property2.currentText()[-2:-1].lower()]
        if SP[0]==SP[1]:
            self._lbl_Warning.setText("Warning:  You cannot specify the same property twice.")
        else:
            self._lbl_Warning.setText("")

        f=[float(self._le_Property1.text()),float(self._le_Property2.text())]
        SI=self._rdo_SI.isChecked()
        self.state1.setState(SP[0],SP[1],f[0],f[1],SI)
        self._lbl_StateProperties.setText(self.makeLabel(self.state1))

#endregion

#region function definitions
def main():
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    sys.exit(app.exec_())
    pass
#end region

#region function calls
if __name__=="__main__":
    main()
#endregion