/*
 *  Dense-CPP
 *
 *     Nils Hamel - nils.hamel@bluewin.ch
 *     Copyright (c) 2016-2019 EPFL, HES-SO Valais
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

    /*! \file   Dense.hpp
     *  \author Nils Hamel <nils.hamel@bluewin.ch>
     *
     *  Dense-CPP - main module
     */

/*
    header - inclusion guard
 */

    # ifndef __SV_DENSE__
    # define __SV_DENSE__

/*
    header - internal includes
 */

    # include "Image.h"
    # include "OpticalFlow.h"

/*
    header - external includes
 */

    # include <iostream>
    # include <fstream>
    # include <opencv2/core/core.hpp>
    # include <opencv2/highgui/highgui.hpp>
    # include <Eigen/Dense>

/*
    header - preprocessor definitions
 */

    /* define mathematical constants */
    # define SV_PI ( 3.14159265358979323846264338327950288419716939937510L )

/*
    header - preprocessor macros
 */

/*
    header - type definition
 */

/*
    header - structures
 */

    typedef struct sv_point_struct {

        double x;
        double y;
        double z;

    } sv_point;

/*
    header - function prototypes
 */

    /* *** */

    sv_point sv_convert_cartesian( long const sv_width, long const sv_height, double sv_phi, double sv_theta );

    /* *** */

    void sv_estimation_load( char const * const sv_estimation_path, Eigen::Matrix3d & sv_r12, Eigen::Vector3d & sv_t12, Eigen::Matrix3d & sv_r23, Eigen::Vector3d & sv_t23 );

    /* *** */

    cv::Mat sv_dense_image_load( char const * const sv_image_path );

    /* *** */

    void sv_dense_image_check( cv::Mat const & sv_image, long const sv_width, long const sv_height, long const sv_depth );

    /* *** */

    int sv_dense_flow( cv::Mat & sv_img_a, cv::Mat & sv_img_b, long const sv_width, long const sv_height, long const sv_depth, DImage & sv_flow_u, DImage & sv_flow_v );

    /* *** */

    void sv_match_compute( long const sv_width, long const sv_height, DImage const & sv_flow_21_u, DImage const & sv_flow_21_v, DImage const & sv_flow_23_u, DImage const & sv_flow_23_v, std::vector < sv_point > & sv_mat_1, std::vector < sv_point > & sv_mat_2, std::vector < sv_point > & sv_mat_3 );

    /*! \brief main function
     *
     *  \param  argc Standard parameter
     *  \param  argv Standard parameter
     *
     *  \return Returns standard exit code
     */

    int main( int argc, char ** argv );

/*
    header - inclusion guard
 */

    # endif

