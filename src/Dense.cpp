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

    Eigen::Vector3d sv_convert_cartesian( long const sv_width, long const sv_height, double sv_phi, double sv_theta ) {

        /* returned structure variable */
        Eigen::Vector3d sv_return;

        /* coordinates re-normalisation */
        sv_phi = ( sv_phi / sv_width ) * 2.0 * SV_PI;

        /* coordinates re-normalisation */
        sv_theta = ( ( sv_theta / ( sv_height - 1 ) ) - 0.5 ) * SV_PI;

        /* compute cartesian coordinates */
        sv_return(0) = cos( sv_theta ) * cos( sv_phi );
        sv_return(1) = cos( sv_theta ) * sin( sv_phi );
        sv_return(2) = sin( sv_theta );

        /* return converted coordinates */
        return( sv_return );

    }

    void sv_convert_to_first_frame( std::vector < Eigen::Vector3d > & sv_mat_1, std::vector < Eigen::Vector3d > & sv_mat_2, std::vector < Eigen::Vector3d > & sv_mat_3, Eigen::Vector3d & sv_cen_1, Eigen::Vector3d & sv_cen_2, Eigen::Vector3d & sv_cen_3, Eigen::Matrix3d const & sv_r12, Eigen::Vector3d const & sv_t12, Eigen::Matrix3d const & sv_r23, Eigen::Vector3d const & sv_t23 ) {

        /* parsing directions vectors */
        for ( long sv_parse( 0 ); sv_parse < sv_mat_1.size(); sv_parse ++ ) {

            /* compute direction - middle camera */
            sv_mat_2[sv_parse] = sv_r12.transpose() * sv_mat_2[sv_parse];

            /* compute direction - last camera */
            sv_mat_3[sv_parse] = sv_r12.transpose() * ( sv_r23.transpose() * sv_mat_3[sv_parse] );

        }

        /* compute center - first camera */
        sv_cen_1 = Eigen::Vector3d::Zero();

        /* compute center - middle camera */
        sv_cen_2 = sv_r12.transpose() * sv_t12;

        /* compute center - last camera */
        sv_cen_3 = sv_cen_2 - ( sv_r12.transpose() * sv_r23.transpose() ) * sv_t23;

    }

/*
    source - optimisation methods
 */


    Eigen::Vector3d sv_dense_optimise_intersect( Eigen::Vector3d const & sv_mat_1, Eigen::Vector3d const & sv_mat_2, Eigen::Vector3d const & sv_mat_3, Eigen::Vector3d  const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3 ) {

        /* intermediate matrix */
        Eigen::Matrix3d sv_w1;
        Eigen::Matrix3d sv_w2;
        Eigen::Matrix3d sv_w3;

        /* intermediate vector */
        Eigen::Vector3d sv_q1;
        Eigen::Vector3d sv_q2;
        Eigen::Vector3d sv_q3;

        /* intermediate computation */
        sv_w1 = Eigen::Matrix3d::Identity() - ( sv_mat_1 * sv_mat_1.transpose() );
        sv_q1 = sv_w1 * sv_cen_1;

        /* intermediate computation */
        sv_w2 = Eigen::Matrix3d::Identity() - ( sv_mat_2 * sv_mat_2.transpose() );
        sv_q2 = sv_w2 * sv_cen_2;

        /* intermediate computation */
        sv_w3 = Eigen::Matrix3d::Identity() - ( sv_mat_3 * sv_mat_3.transpose() );
        sv_q3 = sv_w3 * sv_cen_3;

        /* compute intersection */
        return( ( ( sv_w1 + sv_w2 + sv_w3 ).inverse() ) * ( sv_q1 + sv_q2 + sv_q3 ) );

    }

/*
    source - estimation importation
 */

    void sv_estimation_load( char const * const sv_estimation_path, Eigen::Matrix3d & sv_r12, Eigen::Vector3d & sv_t12, Eigen::Matrix3d & sv_r23, Eigen::Vector3d & sv_t23 ) {

        /* reading matrix variable */
        Eigen::VectorXd sv_read( 24 );

        /* stream variable */
        std::fstream sv_stream;

        /* create stream */
        sv_stream.open( sv_estimation_path, std::ios::in );

        /* check stream */
        if ( sv_stream.is_open() == false ) {

            /* display message */
            std::cerr << "scanvan : error : unable to import estimation" << std::endl;

            /* send message */
            exit( 1 );

        }

        /* import estimation parameter */
        for ( long sv_count( 0 ); sv_count < 24; sv_count ++ ) {

            /* import token */
            sv_stream >> sv_read(sv_count);

        }

        /* delete input stream */
        sv_stream.close();

        /* compose estimation matrix */
        sv_r12(0,0) = sv_read( 0); sv_r12(0,1) = sv_read( 1); sv_r12(0,2) = sv_read( 2);
        sv_r12(1,0) = sv_read( 8); sv_r12(1,1) = sv_read( 9); sv_r12(1,2) = sv_read(10);
        sv_r12(2,0) = sv_read(16); sv_r12(2,1) = sv_read(17); sv_r12(2,2) = sv_read(18);

        /* compose estimation matrix */
        sv_r23(0,0) = sv_read( 4); sv_r23(0,1) = sv_read( 5); sv_r23(0,2) = sv_read( 6);
        sv_r23(1,0) = sv_read(12); sv_r23(1,1) = sv_read(13); sv_r23(1,2) = sv_read(14);
        sv_r23(2,0) = sv_read(20); sv_r23(2,1) = sv_read(21); sv_r23(2,2) = sv_read(22);

        /* compose translation vector */
        sv_t12(0) = sv_read( 3); sv_t12(1) = sv_read(11); sv_t12(2) = sv_read(19);

        /* compose translation vector */
        sv_t23(0) = sv_read( 7); sv_t23(1) = sv_read(15); sv_t23(2) = sv_read(23);

    }

/*
    source - io methods
 */

    void sv_dense_io_export_scene( char const * const sv_path, std::vector < Eigen::Vector3d > const & sv_scene ) {

        /* stream variable */
        std::fstream sv_stream;

        /* create stream */
        sv_stream.open( sv_path, std::ios::out );

        /* check stream */
        if ( sv_stream.is_open() == false ) {

            /* display message */
            std::cerr << "scanvan : error : unable to export scene" << std::endl;

            /* send message */
            exit( 1 );

        }

        /* parsing scene */
        for ( long sv_parse( 0 ); sv_parse < sv_scene.size(); sv_parse ++ ) {

            /* export scene point coordinates */
            sv_stream << sv_scene[sv_parse](0) << " " << sv_scene[sv_parse](1) << " " << sv_scene[sv_parse](2) << std::endl;

        }

        /* close stream */
        sv_stream.close();

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

    void sv_match_compute( long const sv_width, long const sv_height, DImage const & sv_flow_21_u, DImage const & sv_flow_21_v, DImage const & sv_flow_23_u, DImage const & sv_flow_23_v, std::vector < Eigen::Vector3d > & sv_mat_1, std::vector < Eigen::Vector3d > & sv_mat_2, std::vector < Eigen::Vector3d > & sv_mat_3 ) {

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
    source - scene computation
 */

    std::vector < Eigen::Vector3d > sv_dense_scene_compute( std::vector < Eigen::Vector3d > const & sv_mat_1, std::vector < Eigen::Vector3d > const & sv_mat_2, std::vector < Eigen::Vector3d > const & sv_mat_3, Eigen::Vector3d const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3 ) {

        /* returned structure variable */
        std::vector < Eigen::Vector3d > sv_scene;

        /* parsing direction vectors */
        for ( long sv_parse( 0 ); sv_parse < sv_mat_1.size(); sv_parse ++ ) {

            /* compute and push optimised intersection */
            sv_scene.push_back( sv_dense_optimise_intersect( sv_mat_1[sv_parse], sv_mat_2[sv_parse], sv_mat_3[sv_parse], sv_cen_1, sv_cen_2, sv_cen_3 ) );

        }

        /* return computed scene */
        return( sv_scene );

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
        std::vector < Eigen::Vector3d > sv_mat_1;
        std::vector < Eigen::Vector3d > sv_mat_2;
        std::vector < Eigen::Vector3d > sv_mat_3;

        /* center variable */
        Eigen::Vector3d sv_cen_1;
        Eigen::Vector3d sv_cen_2;
        Eigen::Vector3d sv_cen_3;

        /* estimation parameters */
        Eigen::Matrix3d sv_r12;
        Eigen::Matrix3d sv_r23;
        Eigen::Vector3d sv_t12;
        Eigen::Vector3d sv_t23;

        /* scene variable */
        std::vector < Eigen::Vector3d > sv_scene;

        /* check consistency */
        if ( argc != 6 ) {

            /* display message */
            std::cerr << "scanvan : error : wrong usage" << std::endl;

            /* send message */
            return( 1 );

        }

        /* import estimation parameters */
        sv_estimation_load( argv[4], sv_r12, sv_t12, sv_r23, sv_t23 );

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

        /* compute common frame - aligned on first camera */
        sv_convert_to_first_frame( sv_mat_1, sv_mat_2, sv_mat_3, sv_cen_1, sv_cen_2, sv_cen_3, sv_r12, sv_t12, sv_r23, sv_t23 );

        /* compute scene */
        sv_scene = sv_dense_scene_compute( sv_mat_1, sv_mat_2, sv_mat_3, sv_cen_1, sv_cen_2, sv_cen_3 );

        /* export computed scene */
        sv_dense_io_export_scene( argv[5], sv_scene );

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

