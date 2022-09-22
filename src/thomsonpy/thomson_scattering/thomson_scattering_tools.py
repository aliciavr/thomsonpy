"""
.. module:: thomson_scattering_tools
        :platform: Unix
        :synopsis: tools for computing the Thomson scattering across the line of sight 
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>
"""
from math import *
from matplotlib import plt
import numpy as np
import thomsonpy.constants.universal_constants as uc
from thomsonpy.constants.units import *
import thomsonpy.config.thomson_scattering_params as tsp
import thomsonpy.thomson_scattering.ne_models as ne

RSOL = RSOL_TO_METERS


class ThomsonGeometry:
    """
    This class manages the Thomson scattering geometry and also offers many methods useful for frequent calculations
    related with the Thomson scattering geometry. Additionally, its internal structure has been designed to support a
    ray-tracing process in the :math:`z` direction between the observer coordinates of view (:math:`O`) and the scattering
    coordinates (:math:`Q`).
    """

    def __init__(self, sun_center, observer, target, radius):
        """
        tangential_intensity creates a ThomsonGeometry object, which stores and maintains the fundamental
        magnitudes of the Thomson scattering geometry.

        :param sun_center: Center of the Sun (S) in coordinates of the Observer Reference
                            System.
        :type sun_center: numpy.ndarray([float, float, float])
        :param target: Scattering Point (Q) in coordinates of the Observer Reference
                          System.
        :type target: numpy.ndarray([float, float, float])
        :param radius: solar radius. 
        :type radius: float

        :return: an object with a coherent geometry according to the Thomson scattering.
        :rtype: thomsonpy.thomson_scattering.ThomsonGeometry
        """
        # Zero margin value for floating and double comparisons.
        self.__zero = 1E-3

        # Radius of the star.
        self.__radius = radius
        # Position of the Observer.
        self.__observer = observer
        # Position of the Center of the Sun.
        self.__sun_center = sun_center

        # Distance OQ between the Observer and the Scattering Point (z).
        self.__OQ_dist = np.linalg.norm(target - self.__observer)
        # Distance OS between the Observer and the Center of the Sun (x).
        self.__OS_dist = np.linalg.norm(sun_center - self.__observer)
        # Distance SQ between the Center of the Sun and the Scattering Point (d).
        self.__SQ_dist = np.linalg.norm(target - sun_center)

        # Elongation between SQ and OQ distances.
        self.__elongation = np.arcsin(self.__SQ_dist / self.__OQ_dist)

        # Unitary directional vector pointing from the Observer to the Scattering 
        # Point.
        v = target - self.__observer
        norm = np.linalg.norm(v)
        if (norm > (0 + self.__zero)) or (norm < (0 - self.__zero)):
            self.__direction = v / norm
        else:
            self.__direction = np.array([0, 0, 0])

    def __str__(self) -> str:
        """
        A brief description of a ThomsonGeometry object.
        
        :return: a brief description.
        :rtype: string
        """
        return 'Radius = {}\nObserver = {}\nSun Center at {}\nOQ = {}\nOS = {}\nSQ = {}\nelongation = {}º\ndirection ' \
               '= {}'.format(self.__radius, self.__observer, self.__sun_center, self.__OQ_dist, self.__OS_dist,
                             self.__SQ_dist, degrees(self.__elongation), self.__direction)

    def __repr__(self) -> str:
        return '{self.__class__.__name__}({self.__radius, self.__observer, self.__sun_center, self.__OQ_dist, ' \
               'self.__OS_dist, self.__SQ_dist, self.__elongation, self.__direction})'.format(self=self)

    def get_target(self, s):
        """
        Computes the target :math:`\\vec{t}` of Thomson scattering (a.k.a. real-time scattering coordinates) from a
        pre-computed unitary direction :math:`\\vec{d}_u` with the given parameter s.
        
        .. math::
            \\vec{t} = \\vec{o} + s \\vec{d}_u
        
        :param s: parameter of the segment.
        :type s: float
        
        :return: a target :math:`\\vec{t}` for Thomson scattering.
        :rtype: np.ndarray([float, float, float])
        """
        target = self.__observer + s * self.__direction
        return target

    def get_radius(self):
        """
        tangential_intensity gets the radius of the star.

        :return: radius of the Sun.
        :rtype: float
        """
        return self.__radius

    def get_distance_os(self):
        """
        tangential_intensity gets the distance between the Observer (O) and the Center of the Sun (S).
        
        :return: distance OS.
        :rtype: float
        """
        return self.__OS_dist

    def get_distance_oq(self):
        """
        tangential_intensity gets the distance between the Observer (O) and the Scattering Point (Q).

        :return: distance OQ.
        :rtype: float
        """
        return self.__OQ_dist

    def get_distance_sq(self):
        """
        tangential_intensity gets the distance between the Center of the Sun (S) and the Scattering Point (Q).

        :return: distance SQ.
        :rtype: float
        """
        return self.__SQ_dist

    def get_elongation(self):
        """
        tangential_intensity gets the elongation angle, :math:`\epsilon`, between the distance OS and
        the distance OQ.

        :return: elongation (rad).
        :rtype: float
        """
        return self.__elongation

    @staticmethod
    def compute_omega_from_d(d):
        """
        tangential_intensity gets a value for the omega angle between the distance QS and the distance QT as a function
        of the distance d and the radius of the Sun.

        .. math::
            \Omega = \\arcsin{\\frac{RSol}{d}}

        :param d: distance QS (distance d).
        :type d: float

        :return: angle omega (rad).
        :rtype: float
        """
        return asin(tsp.SOLAR_RADIUS / d)

    @staticmethod
    def compute_z_from(x, d, chi):
        """
        tangential_intensity gets a value for the distance OQ (z) as a function of the distance OS (x),
        the distance QS (d) and the scattering angle (S-Q-O).

        .. math::
            z = d \cos{\chi} + \sqrt{d^2(\cos{\chi} - 1)^2 + x^2}
        
        :param x: distance OS.
        :type x: float
        :param d: distance QS.
        :type d: float
        :param chi: scattering angle (rad) between distance OQ and distance QS.
        :type chi: float
        
        :return: distance OQ.
        :rtype: float

        """

        return d * cos(chi) + sqrt(pow(d, 2) * (pow(cos(chi), 2) - 1) + pow(x, 2))

    @staticmethod
    def compute_chi_from(x, z, epsilon):
        """
        tangential_intensity gets a value for the scattering angle, S-Q-O (rad) as a function of
        the distance OS (x), the distance OQ (z) and the elongation angle (:math:`\epsilon`), 
        S-Q-O.
        
        Distance QS (d) is needed for the computations, and therefore obtained from the function
        ``thomsonpy.thomson_scattering.thomson_scattering_tools.ThomsonGeometry.compute_d_from(x, z, epsilon)``
        
        .. math::
            \chi = \\arcsin{\\frac{x \sin{\epsilon}}{d}} \quad where \quad 
            d = \sqrt{x^2 + z^2 - 2 x z \cos{\epsilon}}

        :param x: distance OS.
        :type x: float
        :param z: distance OQ.
        :type z: float
        :param epsilon: elongation angle, Q-O-S (rad).
        :type epsilon: float
        
        :return: scattering angle, S-Q-O (rad).
        :rtype: float
        """

        return np.arcsin(x * np.sin(epsilon) / ThomsonGeometry.compute_d_from(x, z, epsilon))

    @staticmethod
    def compute_d_from(x, z, epsilon):
        """
        tangential_intensity gets a value for the distance QS (d) as a function of the distance OS (x),
        the distance OQ (z) and the elongation angle, :math:`\epsilon` (Q-O-S).
        
        .. math::
            d = \sqrt{x^2 + z^2 - 2 x z \cos{\epsilon}}
            
        :param x: distance OS.
        :type x: float
        :param z: distance OQ.
        :type z: float
        :param epsilon: elongation angle, Q-O-S (rad).
        :type epsilon: float
        
        :return: distance QS (d).
        :rtype: float
        """

        return pow(pow(x, 2) + pow(z, 2) - 2 * x * z * cos(epsilon), 0.5)

    @staticmethod
    def compute_omega_from(x, z, epsilon):
        """
        tangential_intensity gets the angle :math:`\Omega` (T-Q-S) as a function of the distance OS (x), the distance OQ
        (z) and the elongation angle, :math:`\epsilon` (Q-O-S).
        
        Functions ``thomsonpy.thomson_scattering.thomson_scattering_tools.ThomsonGeometry.compute_d_from(x, z,
        epsilon)`` and ``thomsonpy.thomson_scattering.thomson_scattering_tools.ThomsonGeometry.compute_omega_from_d(
        d)`` are used in the computations to obtain d and :math:`\Omega`.
        
        .. math::
            \Omega = \\arcsin{\\frac{RSol}{d}} \quad \quad \quad \quad
            where \quad d = \sqrt{x^2 + z^2 - 2 x z \cos{\epsilon}}
        
        :param x: distance OS.
        :type x: float
        :param z: distance OQ.
        :type z: float
        :param epsilon: elongation angle, Q-O-S (rad).
        :type epsilon: float
        
        :return: angle :math:`\Omega`, T-Q-S (rad).
        :rtype: float
        """

        return ThomsonGeometry.compute_omega_from_d(ThomsonGeometry.compute_d_from(x, z, epsilon))

    @staticmethod
    def compute_z(x, epsilon, phi):
        """
        Gets a value for the distance OQ (z) as a function of the distance OS (x),
        the elongation angle :math:`\epsilon`, and angle :math:`\phi` (O-S-Q).

        .. math::
            z = \\frac{(x \\tan \phi)}{(\sin{\epsilon} + \\tan(\phi) \cos{\epsilon})}
        
        :param x: distance OS.
        :type x: float
        :param epsilon: scattering angle (rad).
        :type epsilon: float
        :param phi: angle O-S-Q (rad).
        :type phi: float
        
        :return: distance OQ.
        :rtype: float
        """

        return (x * tan(phi)) / (sin(epsilon) + tan(phi) * cos(epsilon))


def wien_law(temperature):
    """
    Given a value for temperature it gives the value of the wavelength with 
    the maximum black body emission.
    
    .. math::
        \lambda_{max} = \\frac{2.8978E-3}{T}
    
    :param temperature: temperature of the black body.
    :type temperature: float
    
    :return: the wavelength value with the maximum black body emission.
    :rtype: float
    """

    return 2.8978E-3 / temperature


def planck_law(temperature, wave, is_wavelength=True):
    """
    tangential_intensity computes the Planck's law for the black body radiation.
    
    Two versions of this law are implemented:
    
    * The one based on the wavelength:
        .. math::
            I_\lambda = \\frac{2 h c^2}{\lambda^5}\\frac{1}{e^{\\frac{hc}{\lambda k T}} - 1}
    * The one based on the frequency:
        .. math::
            I_\\nu = \\frac{2 h \\nu^3}{c^2}\\frac{1}{e^{\\frac{h \\nu}{k T}}-1}
            
    Where:
    
    * :math:`h = 6.62607015·10^{-34} J·s` Planck's constant.  
    * :math:`\lambda` is the wavelength of the wave.
    * :math:`\\nu` is the frequency of the wave.
    * :math:`k = 1.38·10^{-23}J·K^{-1}` Boltzmann's constant.
    * :math:`c = 3·10^8m·s^{-1}` is speed of light.
    * :math:`T = 5778 K` is temperature of the black body (in this case the Sun).
        
    
    :param temperature: temperature of the black body, in this case, of the Sun.
    :type temperature: float
    
    :param wave: property of the wave (it can be the wavelength :math:`\lambda` or the frequency :math:`\\nu`) defined
    by the parameter ``is_wavelength``.
    :type wave: float
    
    :param is_wavelength: ``True`` if the parameter ``wave`` represents the wavelength :math:`\lambda` of the wave,
    ``False`` if the parameter ``wave`` represents the frequency :math:`\\nu` of the wave. The default value is
    ``True``.
    :type is_wavelength: boolean
    
    :return: radiation of the black body as a function of the wavelength :math:`I_\lambda` or as a function of the
    frequency :math:`I_\\nu`.
    :rtype: float

    """
    i0 = 0
    if is_wavelength:
        # Planck's Law as function of the wavelength.
        i0 = (2 * uc.h * uc.c ** 2) / (wave ** 5 * (exp((uc.h * uc.c) / (wave * uc.k * temperature)) - 1))
    else:
        # Planck's Law as a function of the frequency.
        i0 = (2 * uc.h * wave ** 3) / (uc.c ** 2 * (np.exp((uc.h * wave) / (uc.k * temperature)) - 1))

    return i0


def __allen_clv(wave, theta, check=0):
    """
    Allen's Astrophysical Quantities, Springer, 2000 
    get u(lambda) from Allen tables. I(theta)/I{(O)
    theta = angle between the Sun's radius vector and the line of sight. In rad.
    wave = wavelength in Angstroms
    forth parameter is u!!!!!
    """

    def coefs(wave):
        u = (- 8.9829751 + 0.0069093916 * wave - 1.8144591e-6 * wave ** 2 + 2.2540875e-10 * wave ** 3 -
             1.3389747e-14 * wave ** 4 + 3.0453572e-19 * wave ** 5)
        v = (+ 9.2891180 - 0.0062212632 * wave + 1.5788029e-6 * wave ** 2 - 1.9359644e-10 * wave ** 3 +
             1.1444469e-14 * wave ** 4 - 2.5994940e-19 * wave ** 5)
        return u, v

    try:
        if check == 1:
            gamma = np.arange(0, 90, 1) * np.pi / 180.
            u, v = coefs(wave)
            I = 1 - (u + v) + (u + v) * np.cos(gamma)
            clv = 1 - u - v + u * np.cos(gamma) + v * np.cos(gamma)
            plt.plot(gamma, clv)
            plt.plot(gamma, I, '--')
    except:
        pass
    u, v = coefs(wave)
    if check == 2:
        return u + v
    return 1 - u - v + u * np.cos(theta) + v * np.cos(theta), u, v, u + v


def coefficient_limb_darkening(wavelength, units=1E+10):
    """
    Adaptation from the private function ``__allen_clv`` to obtain the coefficient
    of limb-darkening :math:`u` needed for the Thomson scattering computations.
    
    The ``__allen_clv`` is a function based at the Allen tables (Astrophysical Quantities, 
    Springer 2000) an implemented in the SPG (IAA - CSIC) repositories.

    :param wavelength: wavelength in I.S. units (meters).
    :type wavelength: float
    :param units: unit conversion factor, by default :math:`1E+10`, the one used to convert from meters to Angstroms.
    Angstroms are the unit used in ``__allev_clv`` function. If the parameter ``wavelength`` is no in I.S. units, a new
    value should be given.
    :type units: float
    
    :return: coefficient of limb-darkening :math:`u`.
    :rtype: float
    """

    wavelength = wavelength * units
    return __allen_clv(wavelength, 1, 2)


def van_de_hulst_a(omega):
    """
    tangential_intensity gets the coefficient A of van de Hulst as a function of the angle :math:`\Omega`.
    
    .. math::
        A = \cos{\Omega} \sin^2{\Omega}
    
    :param omega: angle :math:`\Omega`, T-Q-S (rad).
    :type omega: float

    :return: Coefficient A of van de Hulst.
    :rtype: float
    """

    return cos(omega) * pow(sin(omega), 2)


def van_de_hulst_b(omega):
    """
    tangential_intensity gets the coefficient B of van de Hulst as a function of the angle :math:`\Omega`.
    
    .. math::
        B = -\\frac{1}{8} \left(1 - 3 \sin^2{\Omega} - \\frac{\cos^2{\Omega}}{\sin{\Omega}} (1 + 3 \sin^2{\Omega}) \log{\left(\\frac{1 + \sin{\Omega}}{\cos{\Omega}}\\right)}\\right)
    
    :param omega: angle :math:`\Omega`, T-Q-S (rad).
    :type omega: float

    :return: Coefficient B of van de Hulst.
    :rtype: float
    """

    return -1 / 8 * (1 - 3 * pow(sin(omega), 2) - (pow(cos(omega), 2) / sin(omega)
                                                   * (1 + 3 * pow(sin(omega), 2)) * log((1 + sin(omega)) / cos(omega))))


def van_de_hulst_c(omega):
    """
    tangential_intensity gets the coefficient C of van de Hulst as a function of the angle :math:`\Omega`.
    
    .. math::
        C = \\frac{4}{3} - \cos{\Omega} - \\frac{\cos^3{\Omega}}{3}
    
    :param omega: angle :math:`\Omega`, T-Q-S (rad).
    :type omega: float

    :return: Coefficient C of van de Hulst.
    :rtype: float
    """

    return 4 / 3 - cos(omega) - pow(cos(omega), 3) / 3


def van_de_hulst_d(omega):
    """
    tangential_intensity gets the coefficient D of van de Hulst as a function of the angle :math:`\Omega`.
    
    .. math::
        D = \\frac{1}{8}  \left(5 + \sin^2{\Omega} - \\frac{\cos^2{\Omega}}{\sin{\Omega}} (5 - \sin^2{\Omega}) \log{\left(\\frac{1 + \sin{\Omega}}{\cos{\Omega}}\\right)}\\right)
                  
    :param omega: angle :math:`\Omega`, T-Q-S (rad).
    :type omega: float

    :return: Coefficient D of van de Hulst.
    :rtype: float
    """

    return 1 / 8 * (5 + pow(sin(omega), 2) - (pow(cos(omega), 2) / sin(omega)) *
                    (5 - pow(sin(omega), 2)) * log((1 + sin(omega)) / cos(omega)))


def van_de_hulst(omega, coefficient):
    """
    tangential_intensity computes any of the four coefficients of van de Hulst given a value for the
    angle :math:`\Omega` (T-Q-S).
    
    The name of the coefficients can be either in uppercase or lowercase.
    
    :param omega: angle :math:`\Omega`, T-Q-S (rad).
    :type omega: float
    
    :param coefficient: name of the coefficient of van de Hulst (in uppercase or lowercase).
    :type coefficient: string
    
    :return: the value of the asked coefficient.
    :rtype: float
    """

    coefficient = coefficient.upper()

    # Coefficients of van de Hulst as a function of angle omega (T-Q-S).
    coeff_value = -1
    if coefficient == 'A':
        # Coefficient A.
        coeff_value = van_de_hulst_a(omega)
    elif coefficient == 'B':
        # Coefficient B.
        coeff_value = van_de_hulst_b(omega)
    elif coefficient == 'C':
        # Coefficient C.
        coeff_value = van_de_hulst_c(omega)
    elif coefficient == 'D':
        # Coefficient D.
        coeff_value = van_de_hulst_d(omega)
    else:
        # Not valid value for coefficient of van de Hulst.
        print("[vanDeHulstSunDistance_ERR]::Incorrect coefficient.")
    return coeff_value


def polarized_intensity(i0, sigma_e, z, omega, chi, u):
    """
    tangential_intensity computes the polarized intensity :math:`I_P` as a function of the initial intensity
    :math:`I_0`, given by the Planck's law, the cross-section for perpendicular Thomson scattering :math:`\sigma_e`,
    the distance OQ (z), the angle :math:`\Omega` (T-Q-S), the scattering angle :math:`\chi` (S-Q-O) and the
    coefficient of limb-darkening :math:`u`.
        
    .. math::
        I_P = I_0 \\frac{\pi \sigma_e}{2 z^2} \sin^2{\chi} ((1-u) A + u B)
    
    :param i0: initial intensity, given by the Planck's law.
    :type i0: float
    :param sigma_e: cross-section for perpendicular scattering :math:`\sigma_e`.
    :type sigma_e: float
    :param z: distance OQ.
    :type z: float
    :param omega: angle :math:`\Omega`, T-Q-S. 
    :type omega: float
    :param chi: scattering angle :math:`\chi`, S-Q-O.
    :type chi: float
    :param u: coefficient of limb-darkening.
    :type u: float
    
    :return: polarized intensity :math:`I_P`.
    :rtype: float
    """

    return i0 * (pi * sigma_e) / (2 * pow(z, 2)) * pow(sin(chi), 2) * (
            (1 - u) * van_de_hulst(omega, 'A') + u * van_de_hulst(omega, 'B'))


def tangential_intensity(i0, sigma_e, z, omega, u):
    """
    tangential_intensity computes the tangential intensity :math:`I_T` as a function of the initial intensity
    :math:`I_0`, given by the Planck's law, the cross-section for perpendicular Thomson scattering :math:`\sigma_e`,
    the distance OQ (z), the angle :math:`Omega` (T-Q-S) and the coefficient of limb-darkening :math:`u`.

    .. math::
        I_T = I_0 \\frac{\pi \sigma_e}{2 z^2} ((1 - u) C + u D)

    :param i0: initial intensity, given by the Planck's law.
    :type i0: float
    :param sigma_e: cross-section for perpendicular scattering :math:`\sigma_e`.
    :type sigma_e: float
    :param z: distance OQ.
    :type z: float
    :param omega: angle :math:`\Omega`, T-Q-S.
    :type omega: float
    :param u: coefficient of limb-darkening.
    :type u: float

    :return: tangential intensity :math:`I_T`.
    :rtype: float
    """

    return i0 * (pi * sigma_e) / (2 * pow(z, 2)) * ((1 - u) * van_de_hulst(omega, 'C') + u * van_de_hulst(omega, 'D'))


def radial_intensity(i0, sigma_e, z, omega, chi, u):
    """
    tangential_intensity computes the radial intensity :math:`I_P` as a function of the initial intensity
    :math:`I_0`, given by the Planck's law, the cross-section for perpendicular Thomson scattering :math:`\sigma_e`,
    the distance OQ (z), the angle :math:`\Omega` (T-Q-S), the scattering angle :math:`\chi` (S-Q-O) and the
    coefficient of limb-darkening :math:`u`.
        
    .. math::
        I_R = I_T - I_P
        
    where
    
    .. math::
        I_T = I_0 \\frac{\pi \sigma_e}{2 z^2} ((1 - u) C + u D) \\
        
        I_P = I_0 \\frac{\pi \sigma_e}{2 z^2} \sin^2{\chi} ((1-u) A + u B)
        
        
    :param i0: initial intensity, given by the Planck's law.
    :type i0: float
    :param sigma_e: cross-section for perpendicular scattering :math:`\sigma_e`.
    :type sigma_e: float
    :param z: distance OQ.
    :type z: float
    :param omega: angle :math:`\Omega`, T-Q-S. 
    :type omega: float
    :param chi: scattering angle :math:`\chi`, S-Q-O.
    :type chi: float
    :param u: coefficient of limb-darkening.
    :type u: float
    
    :return: radial intensity :math:`I_R`.
    :rtype: float
    """
    return tangential_intensity(i0, sigma_e, z, omega, u) - polarized_intensity(i0, sigma_e, z, omega, chi, u)


def total_intensity(i0, sigma_e, z, omega, chi, u):
    """
    tangential_intensity computes the total intensity :math:`I_{Total}` as a function of the initial intensity
    :math:`I_0`, given by the Planck's law, the cross-section for perpendicular Thomson scattering :math:`\sigma_e`,
    the distance OQ (z), the angle :math:`\Omega` (T-Q-S), the scattering angle :math:`\chi` (S-Q-O) and the
    coefficient of limb-darkening :math:`u`.

    The functions :py:function:`thomsonpy.thomson_scattering.thomson_scattering_tools.tangential_intensity` and
    :py:function:`thomsonpy.thomson_scattering.thomson_scattering_tools.polarized_intensity` are used for the
    computations.
    
    .. math::
        I_{Total} = 2I_T - I_P

    :param i0: initial intensity, given by the Planck's law.
    :type i0: float
    :param sigma_e: cross-section for perpendicular scattering :math:`\sigma_e`.
    :type sigma_e: float
    :param z: distance OQ.
    :type z: float
    :param omega: angle :math:`\Omega`, T-Q-S. 
    :type omega: float
    :param chi: scattering angle :math:`\chi`, S-Q-O.
    :type chi: float
    :param u: coefficient of limb-darkening.
    :type u: float
    
    :return: total intensity :math:`I_{Total}`.
    :rtype: float
    """
    return 2 * tangential_intensity(i0, sigma_e, z, omega, u) - polarized_intensity(i0, sigma_e, z, omega, chi, u)


def g_tangential(omega, u, z=1):
    """
    tangential_intensity computes the Thomson scattering factor for the tangential intensity.
    
    .. math::
        G_T = \\frac{\pi \sigma_e}{2 z^2} \left((1-u) C + u D \\right)
        
    :param omega: angle :math:`\Omega`, T-Q-S. 
    :type omega: float
    :param u: coefficient of limb-darkening.
    :type u: float
    :param z: distance OQ, by default is z = 1 for a simplified form of the computations.
    :type z: float
    
    :return: tangential Thomson scattering factor.
    :rtype: float
    """

    return (pi * tsp.SIGMA_E) / (2 * z ** 2) * ((1 - u) * van_de_hulst(omega, 'C') + u * van_de_hulst(omega, 'D'))


def g_polarized(omega, chi, u, z=1):
    """
    tangential_intensity computes the Thomson scattering factor for the polarized intensity.
    
    .. math::
        G_P = \sin^2{\chi} \\frac{\pi \sigma_e}{2 z^2} \left((1 - u) A + u D\\right)
        
    :param omega: angle :math:`\Omega`, T-Q-S. 
    :type omega: float
    :param chi: scattering angle :math:`\chi`, S-Q-O.
    :type chi: float
    :param u: coefficient of limb-darkening.
    :type u: float
    :param z: distance OQ, by default is z = 1 for a simplified form of the computations.
    :type z: float
    
    :return: polarized Thomson scattering factor.
    :rtype: float
    """

    return pow(np.sin(chi), 2) * (pi * tsp.SIGMA_E) / (2 * z ** 2) * (
                (1 - u) * van_de_hulst(omega, 'A') + u * van_de_hulst(omega, 'B'))


def g_radial(omega, chi, u, z=1):
    """
    tangential_intensity computes the Thomson scattering factor for the radial intensity.
    
    .. math::
        G_R = G_T - G_P
        
    :param omega: angle :math:`\Omega`, T-Q-S. 
    :type omega: float
    :param chi: scattering angle :math:`\chi`, S-Q-O.
    :type chi: float
    :param u: coefficient of limb-darkening.
    :type u: float
    :param z: distance OQ, by default is z = 1 for a simplified form of the computations.
    :type z: float
    
    :return: radial Thomson scattering factor.
    :rtype: float
    """

    return g_tangential(omega, u, z) - g_polarized(omega, chi, u, z)


def g_total(omega, chi, u, z=1):
    """
    tangential_intensity computes the Thomson scattering factor for the total intensity.
    
    .. math::
        G_{Total} = 2 G_T - G_P
        
    :param omega: angle :math:`\Omega`, T-Q-S. 
    :type omega: float
    :param chi: scattering angle :math:`\chi`, S-Q-O.
    :type chi: float
    :param u: coefficient of limb-darkening.
    :type u: float
    :param z: distance OQ, by default is z = 1 for a simplified form of the computations.
    :type z: float
    
    :return: radial Thomson scattering factor.
    :rtype: float
    """

    return 2 * g_tangential(omega, u, z) - g_polarized(omega, chi, u, z)


def f_irec_z(x, epsilon, z, u, tg=None, ne_model=None):
    """
    tangential_intensity computes the factor of intensity received after the Thomson scattering phenomenon at a
    coordinates :math:`z_i` of line of sight.
    
    .. math::
        f_{scattering_z}(z) = \\rho(z) G_{Total}(z)
    
    :param x: distance OS.
    :type x: float
    :param epsilon: elongation angle :math:`\epsilon`, Q-O-S.
    :type epsilon: float
    :param z: distance OQ (independent variable of the function).
    :type z: float
    :param u: coefficient of limb-darkening.
    :type u: float
    :param tg: object `ThomsonGeometry` storing information about the geometry applied to this function when working
    with models of electron density and ray - tracing computations.
    :type tg: ThomsonGeometry
    :param ne_model: model of electron density :math:`\\rho` needed for some simulations, not set by default.
    :type ne_model: Octree
    
    :return: factor of intensity received from the Thomson scattering phenomenon at the given coordinates.
    :rtype: float
    """

    d = ThomsonGeometry.compute_d_from(x, z, epsilon)
    omega = ThomsonGeometry.compute_omega_from(x, z, epsilon)
    chi = ThomsonGeometry.compute_chi_from(x, z, epsilon)
    scatt_factor = g_total(omega, chi, u)
    # ne_value = ne.crammer_model(d / tsp.SOLAR_RADIUS)
    ne_value = ne.predictive_science_model(z, tg, ne_model)
    return ne_value * scatt_factor


def irec_z(x, epsilon, ini_z, fin_z, incr_z, u, tg=None, ne_model=None):
    """
    Numerical integration of the Thomson scattering along the line of sight over z.
    
    .. math::
        \int_{z=z_0}^{z=z_f} f_{scattering}(z) dz = |f_{scattering}(z_{i+1}) - f_{scattering}(z_{i})| (z_{i+1} - z_i)
    
    :param x: distance OS.
    :type x: float
    :param epsilon: elongation angle :math:`\epsilon`, Q-O-S.
    :type epsilon: float
    :param ini_z: starting value for the independent variable :math:`z` representing the distance OQ.
    :type ini_z: float
    :param fin_z: final value for the independent variable :math:`z` representing the distance OQ.
    :type fin_z: float
    :param incr_z: step value in the numerical integration for the independent variable :math:`z` representing the
    distance OQ.
    :type incr_z: float
    :param u: coefficient of limb-darkening.
    :type u: float
    :param tg: object `ThomsonGeometry` storing information about the geometry applied to this function when working
    with models of electron density and ray - tracing computations.
    :type tg: ThomsonGeometry
    :param ne_model: model of electron density :math:`\\rho` needed for some simulations, not set by default.
    :type ne_model: Octree
    
    :return: integrated Thomson scattering factor por a given geometry over z.
    :rtype: float
    """

    steps_z = np.arange(ini_z, fin_z + incr_z, incr_z)
    value_irec = 0
    for i in steps_z:
        value_irec += abs(
            f_irec_z(x, epsilon, i + 1, u, tg, ne_model) - f_irec_z(x, epsilon, i, u, tg, ne_model)) * incr_z
    return value_irec


def get_scattered_light(wave, temperature, x, epsilon, ini_z, fin_z, incr_z, tg=None, ne_model=None):
    """  
    Computation of the scattered light by the Thomson Scattering of the K-Corona.
    
    .. math::
        I_{scattered} = I_0 \int_{z=z_0}^{z=z_f} f_{scattering}(z) dz
        
    :param wave: wavelength of the solar spectrum.
    :type wave: float
    :param temperature: temperature of the Sun.
    :type temperature:
    :param x: distance OS.
    :type x: float
    :param epsilon: elongation angle :math:`\epsilon`, Q-O-S.
    :type epsilon: float
    :param ini_z: starting value for the independent variable :math:`z` representing the distance OQ.
    :type ini_z: float
    :param fin_z: final value for the independent variable :math:`z` representing the distance OQ.
    :type fin_z: float
    :param incr_z: step value in the numerical integration for the independent variable :math:`z` representing the
    distance OQ.
    :type incr_z: float
    :param tg: object `ThomsonGeometry` storing information about the geometry applied to this function when working
    with models of electron density and ray - tracing computations.
    :type tg: ThomsonGeometry
    :param ne_model: model of electron density :math:`\\rho` needed for some simulations, not set by default.
    :type ne_model: float
    
    :return: scattered light for a given wavelength, temperature and line of sight (:math:`W·m^-3` in I.S. units).
    :rtype: float
    """
    u = coefficient_limb_darkening(wave)
    i0 = planck_law(temperature, wave)
    scattering = irec_z(x, epsilon, ini_z, fin_z, incr_z, u, tg, ne_model)

    scattered_light = i0 * scattering

    return scattered_light
