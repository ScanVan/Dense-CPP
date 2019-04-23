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

    # include "Dense.hpp"

/*
    source - optical flow
 */

    int sv_dense_compute_flow( cv::Mat & sv_img_a, cv::Mat & sv_img_b, long const sv_width, long const sv_height, long const sv_depth ) {

        /* image length variable */
        long sv_length( sv_width * sv_height * sv_depth );

        /* image variable */
        DImage sv_dimg_a;

        /* image variable */
        DImage sv_dimg_b;

        /* flow variable */
        DImage sv_flow_u;

        /* flow variable */
        DImage sv_flow_v;

        /* image variable */
        DImage sv_warp;

        /* allocate image memory */
        sv_dimg_a.allocate( sv_width, sv_height, sv_depth );

        /* assign color type */
        sv_dimg_a.setColorType( 0 );

        /* allocate image memory */
        sv_dimg_b.allocate( sv_width, sv_height, sv_depth );

        /* assign color type */
        sv_dimg_b.setColorType( 0 );

        /* copy image content */
        memcpy( sv_dimg_a.pData, sv_img_a.data, sv_length * sizeof( double ) );

        /* copy image content */
        memcpy( sv_dimg_b.pData, sv_img_b.data, sv_length * sizeof( double ) );

        /* compute optical flow */
        OpticalFlow::Coarse2FineFlow( sv_flow_u, sv_flow_v, sv_warp, sv_dimg_a, sv_dimg_b, 0.012, 0.75, 20, 7, 1, 30 );

        /* release image memory */
        sv_dimg_a.clear();

        /* release image memory */
        sv_dimg_b.clear();

    }

/*
    source - main function
 */

    int main( int argc, char ** argv ) {

        /* matrix variable */
        cv::Mat sv_img_middle;

        /* matrix variable */
        cv::Mat sv_img_prev;

        /* matrix variable */
        cv::Mat sv_img_next;

        /* check consistency */
        if ( argc != 4 ) {

            /* display message */
            std::cerr << "scanvan : error : wrong usage" << std::endl;

            /* send message */
            return( 1 );

        }

        /* import and check image */
        if ( ! ( sv_img_prev = cv::imread( argv[1], CV_LOAD_IMAGE_COLOR ) ).data ) {

            /* display message */
            std::cerr << "scanvan : error : unable to import image" << std::endl;

            /* send message */
            return( 1 );

        }

        /* convert image matrix to double */
        sv_img_prev.convertTo( sv_img_prev, CV_64FC3 );

        /* image renormalisation */
        sv_img_prev /= 255.0;

        /* import and check image */
        if ( ! ( sv_img_middle = cv::imread( argv[2], CV_LOAD_IMAGE_COLOR ) ).data ) {

            /* display message */
            std::cerr << "scanvan : error : unable to import image" << std::endl;

            /* send message */
            return( 1 );

        }

        /* convert image matrix to double */
        sv_img_middle.convertTo( sv_img_middle, CV_64FC3 );

        /* image renormalisation */
        sv_img_middle /= 255.0;

        /* import and check image */
        if ( ! ( sv_img_next = cv::imread( argv[3], CV_LOAD_IMAGE_COLOR ) ).data ) {

            /* display message */
            std::cerr << "scanvan : error : unable to import image" << std::endl;

            /* send message */
            return( 1 );

        }

        /* convert image matrix to double */
        sv_img_next.convertTo( sv_img_next, CV_64FC3 );

        /* image renormalisation */
        sv_img_next /= 255.0;

        // DEBUG // testing //
        sv_dense_compute_flow( sv_img_prev, sv_img_middle, sv_img_prev.cols, sv_img_prev.rows, sv_img_prev.channels() );

        /* send message */
        return( 0 );

    }
