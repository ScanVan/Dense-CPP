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
    source - geometric methods
 */

    Eigen::Vector3d sv_dense_geometry_cartesian( long const sv_width, long const sv_height, double sv_phi, double sv_theta ) {

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

    void sv_dense_geometry_common( std::vector < Eigen::Vector3d > & sv_mat_1, std::vector < Eigen::Vector3d > & sv_mat_2, std::vector < Eigen::Vector3d > & sv_mat_3, std::vector < Eigen::Vector3d > & sv_mat_4, std::vector < Eigen::Vector3d > & sv_mat_5, Eigen::Vector3d & sv_cen_1, Eigen::Vector3d & sv_cen_2, Eigen::Vector3d & sv_cen_3, Eigen::Vector3d & sv_cen_4, Eigen::Vector3d & sv_cen_5, Eigen::Matrix3d const & sv_r01, Eigen::Vector3d const & sv_t01, Eigen::Matrix3d const & sv_r12, Eigen::Vector3d const & sv_t12, Eigen::Matrix3d const & sv_r23, Eigen::Vector3d const & sv_t23, Eigen::Matrix3d const & sv_r34, Eigen::Vector3d const & sv_t34 ) {

        /* parsing directions vectors */
        for ( long sv_parse( 0 ); sv_parse < sv_mat_1.size(); sv_parse ++ ) {

            /* compute direction - second camera */
            sv_mat_2[sv_parse] = sv_r01.transpose() * sv_mat_2[sv_parse];

            /* compute direction - third camera */
            sv_mat_3[sv_parse] = sv_r01.transpose() * ( sv_r12.transpose() * sv_mat_3[sv_parse] );

            /* compute direction - fourth camera */
            sv_mat_4[sv_parse] = sv_r01.transpose() * ( sv_r12.transpose() * ( sv_r23.transpose() * sv_mat_4[sv_parse] ) );

            /* compute direction - fifth camera */
            sv_mat_5[sv_parse] = sv_r01.transpose() * ( sv_r12.transpose() * ( sv_r23.transpose() * ( sv_r34.transpose() * sv_mat_5[sv_parse] ) ) );

        }

        /* compute center - first camera */
        sv_cen_1 = Eigen::Vector3d::Zero();

        /* compute center - second camera */
        sv_cen_2 = - sv_r01.transpose() * sv_t01;

        /* compute center - third camera */
        sv_cen_3 = sv_cen_2 - ( sv_r01.transpose() * sv_r12.transpose() ) * sv_t12;

        /* compute center - fourth camera */
        sv_cen_4 = sv_cen_3 - ( sv_r01.transpose() * sv_r12.transpose() * sv_r23.transpose() ) * sv_t23;

        /* compute center - fourth camera */
        sv_cen_5 = sv_cen_4 - ( sv_r01.transpose() * sv_r12.transpose() * sv_r23.transpose() * sv_r34.transpose() ) * sv_t34;

    }

    Eigen::Vector3d sv_dense_geometry_intersect( Eigen::Vector3d const & sv_mat_1, Eigen::Vector3d const & sv_mat_2, Eigen::Vector3d const & sv_mat_3, Eigen::Vector3d const & sv_mat_4, Eigen::Vector3d const & sv_mat_5, Eigen::Vector3d  const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3, Eigen::Vector3d const & sv_cen_4, Eigen::Vector3d const & sv_cen_5 ) {

        /* intermediate matrix */
        Eigen::Matrix3d sv_w1;
        Eigen::Matrix3d sv_w2;
        Eigen::Matrix3d sv_w3;
        Eigen::Matrix3d sv_w4;
        Eigen::Matrix3d sv_w5;

        /* intermediate vector */
        Eigen::Vector3d sv_q1;
        Eigen::Vector3d sv_q2;
        Eigen::Vector3d sv_q3;
        Eigen::Vector3d sv_q4;
        Eigen::Vector3d sv_q5;

        /* intermediate computation */
        sv_w1 = Eigen::Matrix3d::Identity() - ( sv_mat_1 * sv_mat_1.transpose() );
        sv_q1 = sv_w1 * sv_cen_1;

        /* intermediate computation */
        sv_w2 = Eigen::Matrix3d::Identity() - ( sv_mat_2 * sv_mat_2.transpose() );
        sv_q2 = sv_w2 * sv_cen_2;

        /* intermediate computation */
        sv_w3 = Eigen::Matrix3d::Identity() - ( sv_mat_3 * sv_mat_3.transpose() );
        sv_q3 = sv_w3 * sv_cen_3;

        /* intermediate computation */
        sv_w4 = Eigen::Matrix3d::Identity() - ( sv_mat_4 * sv_mat_4.transpose() );
        sv_q4 = sv_w4 * sv_cen_4;

        /* intermediate computation */
        sv_w5 = Eigen::Matrix3d::Identity() - ( sv_mat_5 * sv_mat_5.transpose() );
        sv_q5 = sv_w5 * sv_cen_5;

        /* compute intersection */
        return( ( ( sv_w1 + sv_w2 + sv_w3 + sv_w4 + sv_w5 ).inverse() ) * ( sv_q1 + sv_q2 + sv_q3 + sv_q4 + sv_q5 ) );

    }

    double sv_dense_geometry_amplitude( Eigen::Vector3d const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3, Eigen::Vector3d const & sv_cen_4, Eigen::Vector3d const & sv_cen_5 ) {

        /* returned value variable */
        double sv_return( 0.0 );

        /* distances variable */
        double sv_d12( ( sv_cen_1 - sv_cen_2 ).norm() );
        double sv_d13( ( sv_cen_1 - sv_cen_3 ).norm() );
        double sv_d14( ( sv_cen_1 - sv_cen_4 ).norm() );
        double sv_d15( ( sv_cen_1 - sv_cen_5 ).norm() );
        double sv_d23( ( sv_cen_2 - sv_cen_3 ).norm() );
        double sv_d24( ( sv_cen_2 - sv_cen_4 ).norm() );
        double sv_d25( ( sv_cen_2 - sv_cen_5 ).norm() );
        double sv_d34( ( sv_cen_3 - sv_cen_4 ).norm() );
        double sv_d35( ( sv_cen_3 - sv_cen_5 ).norm() );
        double sv_d45( ( sv_cen_4 - sv_cen_5 ).norm() );

        /* assume extremum */
        sv_return = sv_d12;

        /* check extremum consistency */
        if ( sv_return < sv_d13 ) sv_return = sv_d13;
        if ( sv_return < sv_d14 ) sv_return = sv_d14;
        if ( sv_return < sv_d15 ) sv_return = sv_d15;
        if ( sv_return < sv_d23 ) sv_return = sv_d23;
        if ( sv_return < sv_d24 ) sv_return = sv_d24;
        if ( sv_return < sv_d25 ) sv_return = sv_d25;
        if ( sv_return < sv_d34 ) sv_return = sv_d34;
        if ( sv_return < sv_d35 ) sv_return = sv_d35;
        if ( sv_return < sv_d45 ) sv_return = sv_d45;

        /* return amplitude */
        return( sv_return );

    }

/*
    source - i/o methods
 */

    void sv_dense_io_pose( char const * const sv_estimation_path, Eigen::Matrix3d & sv_r01, Eigen::Vector3d & sv_t01, Eigen::Matrix3d & sv_r12, Eigen::Vector3d & sv_t12, Eigen::Matrix3d & sv_r23, Eigen::Vector3d & sv_t23, Eigen::Matrix3d & sv_r34, Eigen::Vector3d & sv_t34 ) {

        /* reading matrix variable */
        Eigen::VectorXd sv_read( 48 );

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
        for ( long sv_count( 0 ); sv_count < 48; sv_count ++ ) {

            /* import token */
            sv_stream >> sv_read(sv_count);

        }

        /* delete input stream */
        sv_stream.close();

        /* compose estimation matrix */
        sv_r01(0,0) = sv_read( 0); sv_r01(0,1) = sv_read( 1); sv_r01(0,2) = sv_read( 2);
        sv_r01(1,0) = sv_read(16); sv_r01(1,1) = sv_read(17); sv_r01(1,2) = sv_read(18);
        sv_r01(2,0) = sv_read(32); sv_r01(2,1) = sv_read(33); sv_r01(2,2) = sv_read(34);

        /* compose estimation matrix */
        sv_r12(0,0) = sv_read( 4); sv_r12(0,1) = sv_read( 5); sv_r12(0,2) = sv_read( 6);
        sv_r12(1,0) = sv_read(20); sv_r12(1,1) = sv_read(21); sv_r12(1,2) = sv_read(22);
        sv_r12(2,0) = sv_read(36); sv_r12(2,1) = sv_read(37); sv_r12(2,2) = sv_read(38);

        /* compose estimation matrix */
        sv_r23(0,0) = sv_read( 8); sv_r23(0,1) = sv_read( 9); sv_r23(0,2) = sv_read(10);
        sv_r23(1,0) = sv_read(24); sv_r23(1,1) = sv_read(25); sv_r23(1,2) = sv_read(26);
        sv_r23(2,0) = sv_read(40); sv_r23(2,1) = sv_read(41); sv_r23(2,2) = sv_read(42);

        /* compose estimation matrix */
        sv_r34(0,0) = sv_read(12); sv_r34(0,1) = sv_read(13); sv_r34(0,2) = sv_read(14);
        sv_r34(1,0) = sv_read(28); sv_r34(1,1) = sv_read(29); sv_r34(1,2) = sv_read(30);
        sv_r34(2,0) = sv_read(44); sv_r34(2,1) = sv_read(45); sv_r34(2,2) = sv_read(46);

        /* compose translation vector */
        sv_t01(0) = sv_read( 3); sv_t01(1) = sv_read(19); sv_t01(2) = sv_read(35);

        /* compose translation vector */
        sv_t12(0) = sv_read( 7); sv_t12(1) = sv_read(23); sv_t12(2) = sv_read(39);

        /* compose translation vector */
        sv_t23(0) = sv_read(11); sv_t23(1) = sv_read(27); sv_t23(2) = sv_read(43);

        /* compose translation vector */
        sv_t34(0) = sv_read(15); sv_t34(1) = sv_read(31); sv_t34(2) = sv_read(47);

    }

    void sv_dense_io_scene( char const * const sv_path, std::vector < Eigen::Vector3d > const & sv_scene, std::vector < Eigen::Vector3i > const & sv_color ) {

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
            sv_stream << sv_scene[sv_parse](0) << " " << sv_scene[sv_parse](1) << " " << sv_scene[sv_parse](2) << " ";

            /* export scene point color */
            sv_stream << sv_color[sv_parse](0) << " " << sv_color[sv_parse](1) << " " << sv_color[sv_parse](2) << std::endl;

        }

        /* close stream */
        sv_stream.close();

    }

    cv::Mat sv_dense_io_image( char const * const sv_path, double const sv_scale ) {

        /* matrix variable */
        cv::Mat sv_import;

        /* matrix variable */
        cv::Mat sv_image;

        /* import and check image */
        if ( ! ( sv_import = cv::imread( sv_path, cv::IMREAD_COLOR ) ).data ) {

            /* display message */
            std::cerr << "scanvan : error : unable to import image" << std::endl;

            /* send message */
            exit( 1 );

        }

        /* resize image */
        cv::resize( sv_import, sv_image, cv::Size(), sv_scale, sv_scale, cv::INTER_AREA );

        /* convert image to double */
        sv_image.convertTo( sv_image, CV_64FC3 );

        /* image renormalisation */
        sv_image /= 255.0;

        /* return read image */
        return( sv_image );

    }

    cv::Mat sv_dense_io_mask( char const * const sv_path, double const sv_scale ) {

        /* matrix variable */
        cv::Mat sv_import;

        /* matrix variable */
        cv::Mat sv_mask;

        /* import and check image */
        if ( ! ( sv_import = cv::imread( sv_path, cv::IMREAD_GRAYSCALE ) ).data ) {

            /* display message */
            std::cerr << "scanvan : error : unable to import mask" << std::endl;

            /* send message */
            exit( 1 );

        }

        /* resize image */
        cv::resize( sv_import, sv_mask, cv::Size(), sv_scale, sv_scale, cv::INTER_NEAREST );

        /* return imported image */
        return( sv_mask );

    }

/*
    source - consistency methods
 */

    void sv_dense_consistent_image( cv::Mat const & sv_image, long const sv_width, long const sv_height, long const sv_depth ) {

        /* check image */
        if ( ( sv_width != sv_image.cols ) || ( sv_height != sv_image.rows ) || ( sv_depth != sv_image.channels() ) ) {

            /* display message */
            std::cerr << "scanvan : error : image size are different" << std::endl;

            /* send message */
            exit( 1 );

        }

    }

/*
    source - optical flow methods
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

    }

/*
    source - matching methods
 */

    void sv_dense_match( cv::Mat const & sv_image, cv::Mat const & sv_mask, long const sv_width, long const sv_height, DImage const & sv_flow_20_u, DImage const & sv_flow_20_v, DImage const & sv_flow_21_u, DImage const & sv_flow_21_v, DImage const & sv_flow_23_u, DImage const & sv_flow_23_v, DImage const & sv_flow_24_u, DImage const & sv_flow_24_v, std::vector < Eigen::Vector3d > & sv_mat_1, std::vector < Eigen::Vector3d > & sv_mat_2, std::vector < Eigen::Vector3d > & sv_mat_3, std::vector < Eigen::Vector3d > & sv_mat_4, std::vector < Eigen::Vector3d > & sv_mat_5, std::vector < Eigen::Vector3i > & sv_color ) {

        /* parsing pointer variable */
        double * sv_p_20_u( sv_flow_20_u.pData );
        double * sv_p_20_v( sv_flow_20_v.pData );
        double * sv_p_21_u( sv_flow_21_u.pData );
        double * sv_p_21_v( sv_flow_21_v.pData );
        double * sv_p_23_u( sv_flow_23_u.pData );
        double * sv_p_23_v( sv_flow_23_v.pData );
        double * sv_p_24_u( sv_flow_24_u.pData );
        double * sv_p_24_v( sv_flow_24_v.pData );

        /* reset matches array */
        sv_mat_1.clear();
        sv_mat_2.clear();
        sv_mat_3.clear();
        sv_mat_4.clear();
        sv_mat_5.clear();

        /* color vector variable */
        Eigen::Vector3i sv_pixel;

        /* parsing central image pixels */
        for ( long sv_y( 0 ); sv_y < sv_height; sv_y ++ ) {

            /* parsing central image pixels */
            for( long sv_x( 0 ); sv_x < sv_width; sv_x ++ ) {

                /* check mask value */
                if ( sv_mask.at <uchar> ( sv_y, sv_x ) != 0 ) {

                    /* compute and assign match elements */
                    sv_mat_1.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, sv_x + ( * sv_p_20_u ), sv_y + ( * sv_p_20_v ) ) );

                    /* compute and assign match elements */
                    sv_mat_2.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, sv_x + ( * sv_p_21_u ), sv_y + ( * sv_p_21_v ) ) );

                    /* compute and assign match elements */
                    sv_mat_3.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, sv_x, sv_y ) );

                    /* compute and assign match elements */
                    sv_mat_4.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, sv_x + ( * sv_p_23_u ), sv_y + ( * sv_p_23_v ) ) );

                    /* compute and assign match elements */
                    sv_mat_5.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, sv_x + ( * sv_p_24_u ), sv_y + ( * sv_p_24_v ) ) );

                    /* compose color vector */
                    sv_pixel(0) = sv_image.at <cv::Vec3d> ( sv_y, sv_x )[2] * 255.0;
                    sv_pixel(1) = sv_image.at <cv::Vec3d> ( sv_y, sv_x )[1] * 255.0;
                    sv_pixel(2) = sv_image.at <cv::Vec3d> ( sv_y, sv_x )[0] * 255.0;

                    /* push color vector */
                    sv_color.push_back( sv_pixel );

                }

                /* update pointers */
                sv_p_20_u ++;
                sv_p_20_v ++;
                sv_p_21_u ++;
                sv_p_21_v ++;
                sv_p_23_u ++;
                sv_p_23_v ++;
                sv_p_24_u ++;
                sv_p_24_v ++;

            }

        }


    }

/*
    source - scene methods
 */

    std::vector < Eigen::Vector3d > sv_dense_scene( std::vector < Eigen::Vector3d > const & sv_mat_1, std::vector < Eigen::Vector3d > const & sv_mat_2, std::vector < Eigen::Vector3d > const & sv_mat_3, std::vector < Eigen::Vector3d > const & sv_mat_4, std::vector < Eigen::Vector3d > const & sv_mat_5, Eigen::Vector3d const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3, Eigen::Vector3d const & sv_cen_4, Eigen::Vector3d const & sv_cen_5 ) {

        /* returned structure variable */
        std::vector < Eigen::Vector3d > sv_scene;

        /* parsing direction vectors */
        for ( long sv_parse( 0 ); sv_parse < sv_mat_1.size(); sv_parse ++ ) {

            /* compute and push optimised intersection */
            sv_scene.push_back( sv_dense_geometry_intersect( sv_mat_1[sv_parse], sv_mat_2[sv_parse], sv_mat_3[sv_parse], sv_mat_4[sv_parse], sv_mat_5[sv_parse], sv_cen_1, sv_cen_2, sv_cen_3, sv_cen_4, sv_cen_5 ) );

        }

        /* return computed scene */
        return( sv_scene );

    }

/*
    source - filtering methods
 */

    void sv_dense_filter( double const sv_tol, double const sv_max, std::vector < Eigen::Vector3d > const & sv_scene, std::vector < Eigen::Vector3i > const & sv_color, std::vector < Eigen::Vector3d > & sv_fscene, std::vector < Eigen::Vector3i > & sv_fcolor, std::vector < Eigen::Vector3d > const & sv_mat_1, std::vector < Eigen::Vector3d > const & sv_mat_2, std::vector < Eigen::Vector3d > const & sv_mat_3, std::vector < Eigen::Vector3d > const & sv_mat_4, std::vector < Eigen::Vector3d > const & sv_mat_5, Eigen::Vector3d const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3, Eigen::Vector3d const & sv_cen_4, Eigen::Vector3d const & sv_cen_5 ) {

        /* element radius variable */
        double sv_rad_0( 0.0 );
        double sv_rad_1( 0.0 );
        double sv_rad_2( 0.0 );
        double sv_rad_3( 0.0 );
        double sv_rad_4( 0.0 );
        double sv_rad_5( 0.0 );

        /* element disparity variable */
        double sv_disp_1( 0.0 );
        double sv_disp_2( 0.0 );
        double sv_disp_3( 0.0 );
        double sv_disp_4( 0.0 );
        double sv_disp_5( 0.0 );

        /* compute mean center */
        Eigen::Vector3d sv_cen_0( 0.2 * ( sv_cen_1 + sv_cen_2 + sv_cen_3 + sv_cen_4 + sv_cen_5 ) );

        /* parsing scene elements */
        for ( long sv_parse( 0 ); sv_parse < sv_scene.size(); sv_parse ++ ) {

            /* compute mean point radius */
            sv_rad_0 = ( sv_scene[sv_parse] - sv_cen_0 ).norm();

            /* compute element radius */
            sv_rad_1 = sv_mat_1[sv_parse].transpose() * ( sv_scene[sv_parse] - sv_cen_1 );
            sv_rad_2 = sv_mat_2[sv_parse].transpose() * ( sv_scene[sv_parse] - sv_cen_2 );
            sv_rad_3 = sv_mat_3[sv_parse].transpose() * ( sv_scene[sv_parse] - sv_cen_3 );
            sv_rad_4 = sv_mat_4[sv_parse].transpose() * ( sv_scene[sv_parse] - sv_cen_4 );
            sv_rad_5 = sv_mat_5[sv_parse].transpose() * ( sv_scene[sv_parse] - sv_cen_5 );

            /* compute element disparity */
            sv_disp_1 = ( sv_cen_1 + ( sv_rad_1 * sv_mat_1[sv_parse] ) - sv_scene[sv_parse] ).norm();
            sv_disp_2 = ( sv_cen_2 + ( sv_rad_2 * sv_mat_2[sv_parse] ) - sv_scene[sv_parse] ).norm();
            sv_disp_3 = ( sv_cen_3 + ( sv_rad_3 * sv_mat_3[sv_parse] ) - sv_scene[sv_parse] ).norm();
            sv_disp_4 = ( sv_cen_4 + ( sv_rad_4 * sv_mat_4[sv_parse] ) - sv_scene[sv_parse] ).norm();
            sv_disp_5 = ( sv_cen_5 + ( sv_rad_5 * sv_mat_5[sv_parse] ) - sv_scene[sv_parse] ).norm();

            /* apply filtering condition */
            if ( sv_disp_1 <= sv_tol )
            if ( sv_disp_2 <= sv_tol )
            if ( sv_disp_3 <= sv_tol )
            if ( sv_disp_4 <= sv_tol )
            if ( sv_disp_5 <= sv_tol ) {

                /* apply filtering condition */
                if ( ( sv_rad_0 > 0.0 ) && ( sv_rad_0 < sv_max ) ) {

                    /* element selection - position */
                    sv_fscene.push_back( sv_scene[sv_parse] );

                    /* element selection - color */
                    sv_fcolor.push_back( sv_color[sv_parse] );

                }

            }

        }

    }

/*
    source - main function
 */

    int main( int argc, char ** argv ) {

        /* matrix variable */
        cv::Mat sv_img_pre2;
        cv::Mat sv_img_prev;
        cv::Mat sv_img_middle;
        cv::Mat sv_img_next;
        cv::Mat sv_img_nex2;

        /* matrix variable */
        cv::Mat sv_mask;

        /* reference variable */
        long sv_img_width ( 0 );
        long sv_img_height( 0 );
        long sv_img_depth ( 0 );

        /* flow variable */
        DImage sv_flow_20_u;
        DImage sv_flow_20_v;
        DImage sv_flow_21_u;
        DImage sv_flow_21_v;
        DImage sv_flow_23_u;
        DImage sv_flow_23_v;
        DImage sv_flow_24_u;
        DImage sv_flow_24_v;

        /* matches variable */
        std::vector < Eigen::Vector3d > sv_mat_1;
        std::vector < Eigen::Vector3d > sv_mat_2;
        std::vector < Eigen::Vector3d > sv_mat_3;
        std::vector < Eigen::Vector3d > sv_mat_4;
        std::vector < Eigen::Vector3d > sv_mat_5;

        /* center variable */
        Eigen::Vector3d sv_cen_1;
        Eigen::Vector3d sv_cen_2;
        Eigen::Vector3d sv_cen_3;
        Eigen::Vector3d sv_cen_4;
        Eigen::Vector3d sv_cen_5;

        /* estimation parameters */
        Eigen::Matrix3d sv_r01;
        Eigen::Matrix3d sv_r12;
        Eigen::Matrix3d sv_r23;
        Eigen::Matrix3d sv_r34;
        Eigen::Vector3d sv_t01;
        Eigen::Vector3d sv_t12;
        Eigen::Vector3d sv_t23;
        Eigen::Vector3d sv_t34;

        /* scene variable */
        std::vector < Eigen::Vector3d > sv_scene;
        std::vector < Eigen::Vector3d > sv_fscene;

        /* color variable */
        std::vector < Eigen::Vector3i > sv_color;
        std::vector < Eigen::Vector3i > sv_fcolor;

        /* filtering variable */
        double sv_tol( 0.0 );
        double sv_max( 0.0 );

        /* check consistency */
        if ( argc != 10 ) {

            /* display message */
            std::cerr << "scanvan : error : wrong usage" << std::endl;

            /* display quick help */
            std::cerr << "./Dense [image 1 path] [image 2 path] [image 3 path] [image 4 path] [image 5 path] [pose estimation file] [scene export path] [mask image path] [image scale]" << std::endl;

            /* send message */
            return( 1 );

        }

        /* import estimation parameters */
        sv_dense_io_pose( argv[6], sv_r01, sv_t01, sv_r12, sv_t12, sv_r23, sv_t23, sv_r34, sv_t34 );

        /* import image */
        sv_img_pre2   = sv_dense_io_image( argv[1], atof( argv[9] ) );
        sv_img_prev   = sv_dense_io_image( argv[2], atof( argv[9] ) );
        sv_img_middle = sv_dense_io_image( argv[3], atof( argv[9] ) );
        sv_img_next   = sv_dense_io_image( argv[4], atof( argv[9] ) );
        sv_img_nex2   = sv_dense_io_image( argv[5], atof( argv[9] ) );

        /* extract image reference */
        sv_img_width  = sv_img_middle.cols;
        sv_img_height = sv_img_middle.rows;
        sv_img_depth  = sv_img_middle.channels();

        /* check image */
        sv_dense_consistent_image( sv_img_pre2, sv_img_width, sv_img_height, sv_img_depth );
        sv_dense_consistent_image( sv_img_prev, sv_img_width, sv_img_height, sv_img_depth );
        sv_dense_consistent_image( sv_img_next, sv_img_width, sv_img_height, sv_img_depth );
        sv_dense_consistent_image( sv_img_nex2, sv_img_width, sv_img_height, sv_img_depth );

        /* import mask */
        sv_mask = sv_dense_io_mask( argv[8], atof( argv[9] ) );

    # pragma omp parallel sections
    {

    # pragma omp section
    {

        /* compute optical flow : image 2 -> 0 */
        sv_dense_flow( sv_img_middle, sv_img_pre2, sv_img_width, sv_img_height, sv_img_depth, sv_flow_20_u, sv_flow_20_v );

    }

    # pragma omp section
    {

        /* compute optical flow : image 2 -> 1 */
        sv_dense_flow( sv_img_middle, sv_img_prev, sv_img_width, sv_img_height, sv_img_depth, sv_flow_21_u, sv_flow_21_v );

    }

    } /* parallel sections */

    # pragma omp parallel sections
    {

    # pragma omp section
    {

        /* compute optical flow : image 2 -> 3 */
        sv_dense_flow( sv_img_middle, sv_img_next, sv_img_width, sv_img_height, sv_img_depth, sv_flow_23_u, sv_flow_23_v );

    }

    # pragma omp section
    {

        /* compute optical flow : image 2 -> 4 */
        sv_dense_flow( sv_img_middle, sv_img_nex2, sv_img_width, sv_img_height, sv_img_depth, sv_flow_24_u, sv_flow_24_v );

    }

    } /* parallel sections */

        /* compute matches */
        sv_dense_match( sv_img_middle, sv_mask, sv_img_width, sv_img_height, sv_flow_20_u, sv_flow_20_v, sv_flow_21_u, sv_flow_21_v, sv_flow_23_u, sv_flow_23_v, sv_flow_24_u, sv_flow_24_v, sv_mat_1, sv_mat_2, sv_mat_3, sv_mat_4, sv_mat_5, sv_color );

        /* compute common frame - aligned on first camera */
        sv_dense_geometry_common( sv_mat_1, sv_mat_2, sv_mat_3, sv_mat_4, sv_mat_5, sv_cen_1, sv_cen_2, sv_cen_3, sv_cen_4, sv_cen_5, sv_r01, sv_t01, sv_r12, sv_t12, sv_r23, sv_t23, sv_r34, sv_t34 );

        /* compute scene */
        sv_scene = sv_dense_scene( sv_mat_1, sv_mat_2, sv_mat_3, sv_mat_4, sv_mat_5, sv_cen_1, sv_cen_2, sv_cen_3, sv_cen_4, sv_cen_5 );

        /* compute filtering tolerence values */
        sv_tol = sv_dense_geometry_amplitude( sv_cen_1, sv_cen_2, sv_cen_3, sv_cen_4, sv_cen_5 );
        sv_max = sv_tol *  15.0;
        sv_tol = sv_tol / 150.0;

        /* filter scene */
        sv_dense_filter( sv_tol, sv_max, sv_scene, sv_color, sv_fscene, sv_fcolor, sv_mat_1, sv_mat_2, sv_mat_3, sv_mat_4, sv_mat_5, sv_cen_1, sv_cen_2, sv_cen_3, sv_cen_4, sv_cen_5 );

        /* export computed scene */
        sv_dense_io_scene( argv[7], sv_fscene, sv_fcolor );

        /* display information */
        std::cout << "Remaining points after filtering : " << sv_fscene.size() << " over " << sv_scene.size() << std::endl;

        /* send message */
        return( 0 );

    }

