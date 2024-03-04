import numpy as _np

from splinepy.bezier import Bezier as _Bezier
from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase


class QuadrupelLattice(_TileBase):
    def __init__(self):
        """
        Lattice base cell, consisting of a rectangle with two diagonals in the
        center with all directions having a different thickness associated to them
        """
        self._dim = 2
        self._para_dim = 2
        self._evaluation_points = _np.array(
            [
                [0.5, 0.5],
            ]
        )
        self._n_info_per_eval_point = 4

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        contact_length=0.5,
        **kwargs,  # noqa ARG002
    ):
        """Create a microtile based on the parameters that describe the branch
        thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : np.array
          only first entry is used, defines the internal radii of the
          branches
        parameter_sensitivities: np.ndarray
          correlates with thickness of branches and entouring wall
        contact_length : double
          required for conformity between tiles, sets the length of the center
          block on the tiles boundary

        Returns
        -------
        microtile_list : list(splines)
        """
        min_thickness = 0.05
        max_thickness = 0.1 * (1 + 2**0.5)*2
        if not isinstance(contact_length, float):
            raise ValueError("Invalid Type")
        if not ((contact_length > 0.0) and (contact_length < 1.0)):
            raise ValueError("Contact length must be in (0.,1.)")

        # set to default if nothing is given
        if parameters is None:
            self._logd("Tile request is not parametrized, setting default 0.2")
            parameters = _np.ones((1, 4)) * 0.1
        else:
            if not (
                _np.all(parameters > min_thickness)
                and _np.all(parameters < max_thickness)
            ):
                raise ValueError(
                    "The parameter must be in 0.01 and 1/(1 + sqrt(2))"
                )
            pass
        self.check_params(parameters)

        # Check if user requests derivative splines
        if self.check_param_derivatives(parameter_sensitivities):
            n_derivatives = parameter_sensitivities.shape[2]
            derivatives = []
        else:
            n_derivatives = 0
            derivatives = None

        # Coefficients
        coeff_sqrt_2 = 2**0.5
        coeff_inv_sqrt_2 = 2 ** (-0.5)
        coeff_one_half = 0.5

        splines = []
        for i_derivative in range(n_derivatives + 1):
            # Constant auxiliary values
            if i_derivative == 0:
                cl = contact_length
                thickness_horizontal = coeff_one_half * parameters[0, 0]
                thickness_vertical = coeff_one_half * parameters[0, 1]
                thickness_diag_up = coeff_one_half * parameters[0, 2]
                thickness_diag_down = coeff_one_half * parameters[0, 3]
                v_one_half = 0.5
                v_one = 1.0
                v_zero = 0.0
                alpha = _np.arctan2((v_one_half-thickness_horizontal),
                                    (v_one_half-thickness_vertical))
                diag_y_up = thickness_diag_up/(2*_np.cos(alpha))
                diag_x_up = thickness_diag_up/(2*_np.sin(alpha))
                diag_y_down = thickness_diag_down/(2*_np.cos(alpha))
                diag_x_down = thickness_diag_down/(2*_np.sin(alpha))
            else:
                cl = 0.0
                thickness_horizontal = (
                    coeff_one_half
                    * parameter_sensitivities[0, 0, i_derivative - 1]
                )
                thickness_vertical = (
                    coeff_one_half * parameters[0, 1, i_derivative - 1]
                )
                thickness_diag_up = (
                    coeff_one_half * parameters[0, 2, i_derivative - 1]
                )
                thickness_diag_down = (
                    coeff_one_half * parameters[0, 3, i_derivative - 1]
                )
                v_one_half = 0.0
                v_one = 0.0
                v_zero = 0.0

            # Set variables (with notation x01 = 1 - x11 )
            # Variables in Vertical direction
            v01 = thickness_horizontal
            v02 = thickness_horizontal + diag_y_up
            # v02 = thickness_vertical + coeff_sqrt_2 * thickness_diag_up
            v03 = thickness_horizontal + diag_y_down
            v11 = v_one - v01
            v12 = v_one - v02
            v13 = v_one - v03

            # Variables in Horizontal direction
            h01 = thickness_vertical
            h02 = thickness_vertical + diag_x_up

            # h02 = thickness_horizontal + coeff_sqrt_2 * thickness_diag_up
            h03 = thickness_vertical + diag_x_down
            h11 = v_one - h01
            h12 = v_one - h02
            h13 = v_one - h03

            # Omnidirectional Variables
            c01 = v_zero
            c02 = (v_one - cl) * 0.5

            c03_x = 2.0*thickness_diag_down*_np.sqrt(thickness_horizontal**2 - thickness_horizontal + thickness_vertical**2 - thickness_vertical + 0.5)/(8.0*thickness_horizontal - 4.0) - 2.0*thickness_diag_up*_np.sqrt(thickness_horizontal**2 - thickness_horizontal + thickness_vertical**2 - thickness_vertical + 0.5)/(8.0*thickness_horizontal - 4.0) + 4.0*thickness_horizontal/(8.0*thickness_horizontal - 4.0) - 2.0/(8.0*thickness_horizontal - 4.0)

            c03_y = 2.0*thickness_diag_down*_np.sqrt(thickness_horizontal**2 - thickness_horizontal + thickness_vertical**2 - thickness_vertical + 0.5)/(8.0*thickness_vertical - 4.0) + 2.0*thickness_diag_up*_np.sqrt(thickness_horizontal**2 - thickness_horizontal + thickness_vertical**2 - thickness_vertical + 0.5)/(8.0*thickness_vertical - 4.0) + 4.0*thickness_vertical/(8.0*thickness_vertical - 4.0) - 2.0/(8.0*thickness_vertical - 4.0)


            c04_x = 2.0*thickness_diag_down*_np.sqrt(thickness_horizontal**2 - thickness_horizontal + thickness_vertical**2 - thickness_vertical + 0.5)/(8.0*thickness_horizontal - 4.0) + 2.0*thickness_diag_up*_np.sqrt(thickness_horizontal**2 - thickness_horizontal + thickness_vertical**2 - thickness_vertical + 0.5)/(8.0*thickness_horizontal - 4.0) + 4.0*thickness_horizontal/(8.0*thickness_horizontal - 4.0) - 2.0/(8.0*thickness_horizontal - 4.0)


            c04_y = 2.0*thickness_diag_down*_np.sqrt(thickness_horizontal**2 - thickness_horizontal + thickness_vertical**2 - thickness_vertical + 0.5)/(8.0*thickness_vertical - 4.0) - 2.0*thickness_diag_up*_np.sqrt(thickness_horizontal**2 - thickness_horizontal + thickness_vertical**2 - thickness_vertical + 0.5)/(8.0*thickness_vertical - 4.0) + 4.0*thickness_vertical/(8.0*thickness_vertical - 4.0) - 2.0/(8.0*thickness_vertical - 4.0)
 


            # c05 = (
            #     0.5
            #     * coeff_one_half
            #     * (thickness_horizontal + thickness_vertical)
            # )
            c05_x = thickness_vertical
            c05_y = thickness_horizontal
            c06 = v_one_half
            c15_x = v_one - c05_x
            c15_y = v_one - c05_y
            c14_x = v_one - c04_x
            c14_y = v_one - c04_y
            c13_x = v_one - c03_x
            c13_y = v_one - c03_y
            c12 = v_one - c02
            c11 = v_one

            # Init return value
            spline_list = []

            # 0
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c01, c01],
                        [c05_x, c05_y],
                        [c01, c02],
                        [h01, v02],
                    ],
                )
            )

            # 1
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c01, c01],
                        [c02, c01],
                        [c05_x, c05_y],
                        [h02, v01],
                    ],
                )
            )

            # 2
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c02, c01],
                        [c12, c01],
                        [h02, v01],
                        [h13, v01],
                    ],
                )
            )

            # 3
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c12, c01],
                        [c11, c01],
                        [h13, v01],
                        [c15_x, c05_y],
                    ],
                ),
            )

            # 4
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c15_x, c05_y],
                        [c11, c01],
                        [h11, v03],
                        [c11, c02],
                    ],
                )
            )

            # 5
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h11, v03],
                        [c11, c02],
                        [h11, v12],
                        [c11, c12],
                    ],
                )
            )

            # 6
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h11, v12],
                        [c11, c12],
                        [c15_x, c15_y],
                        [c11, c11],
                    ],
                )
            )
            # 7
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h12, v11],
                        [c15_x, c15_y],
                        [c12, c11],
                        [c11, c11],
                    ],
                )
            )

            # 8
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h03, v11],
                        [h12, v11],
                        [c02, c11],
                        [c12, c11],
                    ],
                )
            )

            # 9
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c05_x, c15_y],
                        [h03, v11],
                        [c01, c11],
                        [c02, c11],
                    ],
                )
            )

            # 10
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c01, c12],
                        [h01, v13],
                        [c01, c11],
                        [c05_x, c15_y],
                    ],
                )
            )

            # 11
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c01, c02],
                        [h01, v02],
                        [c01, c12],
                        [h01, v13],
                    ],
                )
            )

            # 12
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h01, v02],
                        [c05_x, c05_y],
                        [c04_x, c04_y],
                        [c06, c06],
                    ],
                )
            )

            # 13
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c05_x, c05_y],
                        [h02, v01],
                        [c06, c06],
                        [c03_x, c03_y],
                    ],
                )
            )
            
            # 14
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c03_x, c03_y],
                        [h13, v01],
                        [c06, c06],
                        [c15_x, c05_y],
                    ],
                )
            )

            # 15
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c06, c06],
                        [c15_x, c05_y],
                        [c14_x, c14_y],
                        [h11, v03],
                    ],
                )
            )

            # 16
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c06, c06],
                        [c14_x, c14_y],
                        [c15_x, c15_y],
                        [h11, v12],
                    ],
                )
            )

            # 17
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c06, c06],
                        [c15_x, c15_y],
                        [c13_x, c13_y],
                        [h12, v11],
                    ],
                )
            )

            # 18
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [c05_x, c15_y],
                        [c06, c06],
                        [h03, v11],
                        [c13_x, c13_y],
                    ],
                )
            )

            # 19
            spline_list.append(
                _Bezier(
                    degrees=[1, 1],
                    control_points=[
                        [h01, v13],
                        [c04_x, c04_y],
                        [c05_x, c15_y],
                        [c06, c06],
                    ],
                )
            )

            # Pass to output
            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        # Return results
        return (splines, derivatives)
