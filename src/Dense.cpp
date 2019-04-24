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
    source - conversion operation
 */

    sv_point sv_convert_cartesian( long const sv_width, long const sv_height, double sv_phi, double sv_theta ) {

        /* returned structure variable */
        sv_point sv_return;

        /* coordinates re-normalisation */
        sv_phi = ( ( sv_phi - 1 ) / sv_width ) * 2.0 * SV_PI;

        /* coordinates re-normalisation */
        sv_theta = ( ( sv_theta / sv_height ) - 0.5 ) * SV_PI;

        /* compute cartesian coordinates */
        sv_return.x = cos( sv_theta ) * cos( sv_phi );
        sv_return.y = cos( sv_theta ) * sin( sv_phi );
        sv_return.z = sin( sv_theta );

        /* return converted coordinates */
        return( sv_return );

    }

/*
    source - image manipulation
 */

    cv::Mat sv_dense_image_load( char const * const sv_image_path ) {

        /* matrix variable */
        cv::Mat sv_image;

        /* import and check image */
        if ( ! ( sv_image = cv::imread( sv_image_path, CV_LOAD_IMAGE_COLOR ) ).data ) {

            /* display message */
            std::cerr << "scanvan : error : unable to import image" << std::endl;

            /* send message */
            exit( 1 );

        }

        /* convert image to double */
        sv_image.convertTo( sv_image, CV_64FC3 );

        /* image renormalisation */
        sv_image /= 255.0;

        /* return read image */
        return( sv_image );

    }

    void sv_dense_image_check( cv::Mat const & sv_image, long const sv_width, long const sv_height, long const sv_depth ) {

        /* check image */
        if ( ( sv_width != sv_image.cols ) || ( sv_height != sv_image.rows ) || ( sv_depth != sv_image.channels() ) ) {

            /* display message */
            std::cerr << "scanvan : error : image size are different" << std::endl;

            /* send message */
            exit( 1 );

        }

    }

/*
    source - optical flow
 */

    int sv_dense_flow( cv::Mat & sv_img_a, cv::Mat & sv_img_b, long const sv_width, long const sv_height, long const sv_depth, DImage & sv_flow_u, DImage & sv_flow_v ) {

        /* image length variable */
        long sv_length( sv_width * sv_height * sv_depth );

        /* image variable */
        DImage sv_dimg_a;
        DImage sv_dimg_b;

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

        /* release image memory */
        sv_warp.clear();

    }

/*
    source - matches
 */

    void sv_match_compute( long const sv_width, long const sv_height, DImage const & sv_flow_21_u, DImage const & sv_flow_21_v, DImage const & sv_flow_23_u, DImage const & sv_flow_23_v, std::vector < sv_point > & sv_mat_1, std::vector < sv_point > & sv_mat_2, std::vector < sv_point > & sv_mat_3 ) {

        /* parsing pointer variable */
        double * sv_p_21_u( sv_flow_21_u.pData );
        double * sv_p_21_v( sv_flow_21_v.pData );
        double * sv_p_23_u( sv_flow_23_u.pData );
        double * sv_p_23_v( sv_flow_23_v.pData );

        /* reset matches array */
        sv_mat_1.clear();
        sv_mat_2.clear();
        sv_mat_3.clear();

        /* parsing central image pixels */
        for ( long sv_y( 0 ); sv_y < sv_height; sv_y ++ ) {

            /* parsing central image pixels */
            for( long sv_x( 0 ); sv_x < sv_width; sv_x ++ ) {

                /* compute and assign match elements */
                sv_mat_1.push_back( sv_convert_cartesian( sv_width, sv_height, sv_x + ( * sv_p_21_u ), sv_y + ( * sv_p_21_v ) ) );

                /* compute and assign match elements */
                sv_mat_2.push_back( sv_convert_cartesian( sv_width, sv_height, sv_x, sv_y ) );

                /* compute and assign match elements */
                sv_mat_3.push_back( sv_convert_cartesian( sv_width, sv_height, sv_x + ( * sv_p_23_u ), sv_y + ( * sv_p_23_v ) ) );

                /* update pointers */
                sv_p_21_u ++;
                sv_p_21_v ++;
                sv_p_23_u ++;
                sv_p_23_v ++;

            }

        }


    }

/*
    source - main function
 */

    int main( int argc, char ** argv ) {

        /* matrix variable */
        cv::Mat sv_img_prev;
        cv::Mat sv_img_middle;
        cv::Mat sv_img_next;

        /* reference variable */
        long sv_img_width ( 0 );
        long sv_img_height( 0 );
        long sv_img_depth ( 0 );

        /* flow variable */
        DImage sv_flow_21_u;
        DImage sv_flow_21_v;
        DImage sv_flow_23_u;
        DImage sv_flow_23_v;

        /* matches variable */
        std::vector < sv_point > sv_mat_1;
        std::vector < sv_point > sv_mat_2;
        std::vector < sv_point > sv_mat_3;


        /* check consistency */
        if ( argc != 4 ) {

            /* display message */
            std::cerr << "scanvan : error : wrong usage" << std::endl;

            /* send message */
            return( 1 );

        }

        /* import image */
        sv_img_prev   = sv_dense_image_load( argv[1] );
        sv_img_middle = sv_dense_image_load( argv[2] );
        sv_img_next   = sv_dense_image_load( argv[3] );

        /* extract image reference */
        sv_img_width  = sv_img_middle.cols;
        sv_img_height = sv_img_middle.rows;
        sv_img_depth  = sv_img_middle.channels();

        /* check image */
        sv_dense_image_check( sv_img_prev, sv_img_width, sv_img_height, sv_img_depth );
        sv_dense_image_check( sv_img_next, sv_img_width, sv_img_height, sv_img_depth );

        /* compute optical flows : image 2 -> 1 */
        sv_dense_flow( sv_img_middle, sv_img_prev, sv_img_width, sv_img_height, sv_img_depth, sv_flow_21_u, sv_flow_21_v );

        /* compute optical flows : image 2 -> 3 */
        sv_dense_flow( sv_img_middle, sv_img_next, sv_img_width, sv_img_height, sv_img_depth, sv_flow_23_u, sv_flow_23_v );

        /* compute matches */
        sv_match_compute( sv_img_width, sv_img_height, sv_flow_21_u, sv_flow_21_v, sv_flow_23_u, sv_flow_23_v, sv_mat_1, sv_mat_2, sv_mat_3 );

    // DEBUG // CHECK //
# ifdef __DEBUG
    std::fstream __stream; double * __p = sv_flow_u.pData;
    __stream.open( "export_u.dat", std::ios::out );
    for ( int __y( 0 ); __y < sv_img_height; __y ++ ) {
    for ( int __x( 0 ); __x < sv_img_width ; __x ++ ) {
        __stream << * __p << " "; ++__p;
    }
    __stream << std::endl;
    }
    __stream.close();
    // DEBUG // CHECK //
# endif

        /* release flow memory */
        sv_flow_21_u.clear();
        sv_flow_21_v.clear();
        sv_flow_23_u.clear();
        sv_flow_23_v.clear();

        /* send message */
        return( 0 );

    }

