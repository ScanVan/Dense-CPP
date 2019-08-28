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

    void sv_dense_geometry_common( std::vector < Eigen::Vector3d > & sv_mat_1, std::vector < Eigen::Vector3d > & sv_mat_2, std::vector < Eigen::Vector3d > & sv_mat_3, Eigen::Vector3d & sv_cen_1, Eigen::Vector3d & sv_cen_2, Eigen::Vector3d & sv_cen_3, Eigen::Matrix3d const & sv_r12, Eigen::Vector3d const & sv_t12, Eigen::Matrix3d const & sv_r23, Eigen::Vector3d const & sv_t23 ) {

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
        sv_cen_2 = - sv_r12.transpose() * sv_t12;

        /* compute center - last camera */
        sv_cen_3 = sv_cen_2 - ( sv_r12.transpose() * sv_r23.transpose() ) * sv_t23;

    }

    Eigen::Vector3d sv_dense_geometry_intersect( Eigen::Vector3d const & sv_mat_1, Eigen::Vector3d const & sv_mat_2, Eigen::Vector3d const & sv_mat_3, Eigen::Vector3d  const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3 ) {

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

    double sv_dense_geometry_amplitude( Eigen::Vector3d const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3 ) {

        /* returned value variable */
        double sv_return( 0.0 );

        /* distances variable */
        double sv_d12( ( sv_cen_1 - sv_cen_2 ).norm() );
        double sv_d23( ( sv_cen_2 - sv_cen_3 ).norm() );
        double sv_d31( ( sv_cen_3 - sv_cen_1 ).norm() );

        /* assume extremum */
        sv_return = sv_d12;

        /* check extremum consistency */
        if ( sv_return < sv_d23 ) sv_return = sv_d23;
        if ( sv_return < sv_d31 ) sv_return = sv_d31;

        /* return amplitude */
        return( sv_return );

    }

/*
    source - i/o methods
 */

    void sv_dense_io_pose( char const * const sv_estimation_path, Eigen::Matrix3d & sv_r12, Eigen::Vector3d & sv_t12, Eigen::Matrix3d & sv_r23, Eigen::Vector3d & sv_t23 ) {

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

    string type2str(int type) {
      string r;

      uchar depth = type & CV_MAT_DEPTH_MASK;
      uchar chans = 1 + (type >> CV_CN_SHIFT);

      switch ( depth ) {
        case CV_8U:  r = "8U"; break;
        case CV_8S:  r = "8S"; break;
        case CV_16U: r = "16U"; break;
        case CV_16S: r = "16S"; break;
        case CV_32S: r = "32S"; break;
        case CV_32F: r = "32F"; break;
        case CV_64F: r = "64F"; break;
        default:     r = "User"; break;
      }

      r += "C";
      r += (chans+'0');

      return r;
    }

    void sv_dense_flow_nvidia(std::string sv_path, double sv_scale, DImage & sv_flow_u, DImage & sv_flow_v ){
        /* matrix variable */
        cv::Mat sv_import;

        /* matrix variable */
        cv::Mat sv_image;

        /* import and check image */
        if ( ! ( sv_import = cv::imread( sv_path, cv::IMREAD_ANYCOLOR | cv::IMREAD_ANYDEPTH ) ).data) {

            /* display message */
            std::cerr << "scanvan : error : unable to import image" << std::endl;

            /* send message */
            exit( 1 );

        }


        /* resize image */
////        cv::GaussianBlur( sv_import, sv_import, cv::Size(149, 149 ), 0, 0 );
////        cv::namedWindow( "Display window", cv::WINDOW_FREERATIO | cv::WINDOW_GUI_EXPANDED );// Create a window for display.
////        cv::imshow( "Display window", sv_import );
////        cv::waitKey(0);
//        cv::resize( sv_import, sv_image, cv::Size(), sv_scale, sv_scale, cv::INTER_AREA );
//
//
//		cv::GaussianBlur( sv_image, sv_image, cv::Size(31, 31 ), 0, 0 );
//		cv::GaussianBlur( sv_image, sv_image, cv::Size(31, 31 ), 0, 0 );
//
//        sv_flow_u.allocate(sv_image.cols, sv_image.rows, 1);
//        sv_flow_v.allocate(sv_image.cols, sv_image.rows, 1);
//        double * sv_p_u( sv_flow_u.pData );
//        double * sv_p_v( sv_flow_v.pData );
//        for(uint32_t y = 0; y < sv_image.rows;y++){
//            for(uint32_t x = 0; x < sv_image.cols;x++){
//				auto bgr = sv_image.at<cv::Vec3w>(y, x);
//				*sv_p_u = (bgr[2]-32768)/64.0*sv_scale;
//				*sv_p_v = (bgr[1]-32768)/64.0*sv_scale;
//				sv_p_u++;
//				sv_p_v++;
//            }
//        }

        sv_import.convertTo(sv_image, CV_32F, 1.0/64.0, -512.0);

//        cv::namedWindow( "Display window", cv::WINDOW_FREERATIO | cv::WINDOW_GUI_EXPANDED );// Create a window for display.
//        cv::imshow( "Display window", sv_import );

//		cv::GaussianBlur( sv_image, sv_image, cv::Size(31, 31 ), 0, 0 );
//		cv::GaussianBlur( sv_image, sv_image, cv::Size(31, 31 ), 0, 0 );

//		cv::Mat view;
//		sv_image.convertTo(view, CV_16U, 64, 32768);

//        cv::namedWindow( "filtred", cv::WINDOW_FREERATIO | cv::WINDOW_GUI_EXPANDED );// Create a window for display.
//        cv::imshow( "filtred", view );

//        cv::waitKey(0);
		cv::resize( sv_image, sv_image, cv::Size(), sv_scale, sv_scale, cv::INTER_AREA );
		sv_flow_u.allocate(sv_image.cols, sv_image.rows, 1);
		sv_flow_v.allocate(sv_image.cols, sv_image.rows, 1);
		double * sv_p_u( sv_flow_u.pData );
		double * sv_p_v( sv_flow_v.pData );
		for(uint32_t y = 0; y < sv_image.rows;y++){
			for(uint32_t x = 0; x < sv_image.cols;x++){
				auto bgr = sv_image.at<cv::Vec3f>(y, x);
				*sv_p_u = (bgr[2])*sv_scale;
				*sv_p_v = (bgr[1])*sv_scale;
				sv_p_u++;
				sv_p_v++;
			}
		}

    }

/*
    source - matching methods
 */
#ifdef OPTICAL_FLOW_CPP
    void sv_dense_match( cv::Mat const & sv_image, cv::Mat const & sv_mask, long const sv_width, long const sv_height, DImage const & sv_flow_21_u, DImage const & sv_flow_21_v, DImage const & sv_flow_23_u, DImage const & sv_flow_23_v, std::vector < Eigen::Vector3d > & sv_mat_1, std::vector < Eigen::Vector3d > & sv_mat_2, std::vector < Eigen::Vector3d > & sv_mat_3, std::vector < Eigen::Vector3i > & sv_color ) {

        /* parsing pointer variable */
        double * sv_p_21_u( sv_flow_21_u.pData );
        double * sv_p_21_v( sv_flow_21_v.pData );
        double * sv_p_23_u( sv_flow_23_u.pData );
        double * sv_p_23_v( sv_flow_23_v.pData );

        /* reset matches array */
        sv_mat_1.clear();
        sv_mat_2.clear();
        sv_mat_3.clear();

        /* color vector variable */
        Eigen::Vector3i sv_pixel;

        /* parsing central image pixels */
        for ( long sv_y( 0 ); sv_y < sv_height; sv_y ++ ) {

            /* parsing central image pixels */
            for( long sv_x( 0 ); sv_x < sv_width; sv_x ++ ) {

                /* check mask value */
                if ( sv_mask.at <uchar> ( sv_y, sv_x ) != 0 ) {

                    /* compute and assign match elements */
                    sv_mat_1.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, sv_x + ( * sv_p_21_u ), sv_y + ( * sv_p_21_v ) ) );

                    /* compute and assign match elements */
                    sv_mat_2.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, sv_x, sv_y ) );

                    /* compute and assign match elements */
                    sv_mat_3.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, sv_x + ( * sv_p_23_u ), sv_y + ( * sv_p_23_v ) ) );

                    /* compose color vector */
                    sv_pixel(0) = sv_image.at <cv::Vec3d> ( sv_y, sv_x )[2] * 255.0;
                    sv_pixel(1) = sv_image.at <cv::Vec3d> ( sv_y, sv_x )[1] * 255.0;
                    sv_pixel(2) = sv_image.at <cv::Vec3d> ( sv_y, sv_x )[0] * 255.0;

                    /* push color vector */
                    sv_color.push_back( sv_pixel );

                }

                /* update pointers */
                sv_p_21_u ++;
                sv_p_21_v ++;
                sv_p_23_u ++;
                sv_p_23_v ++;

            }

        }


    }
#endif

#ifdef OPTICAL_FLOW_NVIDIA
    double sv_bilinear_sample(double *p, double x, double y, int width){
    	int ix = x;
    	int iy = y;

    	int i00 = iy*width + ix;
    	int i01 = i00 + 1;
    	int i10 = i00 + width;
    	int i11 = i00 + width + 1;

    	double fx = x-ix;
    	double fy = y-iy;

    	return  (p[i00]*(1.0-fx) + p[i01]*fx)*(1.0-fy) + (p[i10]*(1.0-fx) + p[i11]*fx)*fy;
    }

    void sv_dense_match( cv::Mat const & sv_image, cv::Mat const & sv_mask, long const sv_width, long const sv_height, DImage const & sv_flow_12_u, DImage const & sv_flow_12_v, DImage const & sv_flow_23_u, DImage const & sv_flow_23_v, std::vector < Eigen::Vector3d > & sv_mat_1, std::vector < Eigen::Vector3d > & sv_mat_2, std::vector < Eigen::Vector3d > & sv_mat_3, std::vector < Eigen::Vector3i > & sv_color ) {
        /* reset matches array */
        sv_mat_1.clear();
        sv_mat_2.clear();
        sv_mat_3.clear();

        /* color vector variable */
        Eigen::Vector3i sv_pixel;

        /* parsing central image pixels */
        for ( long sv_y( 0 ); sv_y < sv_height; sv_y ++ ) {

            /* parsing central image pixels */
            for( long sv_x( 0 ); sv_x < sv_width; sv_x ++ ) {

                /* check mask value */
                if ( sv_mask.at <uchar> ( sv_y, sv_x ) != 0 ) {
                    double pix_x = sv_x, pix_y = sv_y;

                    /* First picture match */
                    sv_mat_1.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, pix_x, pix_y ) );

                    /* Second picture match */
                    {
						int offset = pix_y*sv_width + pix_x;
						pix_x += sv_flow_12_u.pData[offset];
						pix_y += sv_flow_12_v.pData[offset];
						pix_x = max(0.0, min((double)sv_width-1, pix_x));
						pix_y = max(0.0, min((double)sv_height-1, pix_y));
						sv_mat_2.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, pix_x, pix_y ) );

						/* compose color vector TODO pixel sample rounding*/
	                    sv_pixel(0) = sv_image.at <cv::Vec3d> ( pix_y, pix_x )[2] * 255.0;
	                    sv_pixel(1) = sv_image.at <cv::Vec3d> ( pix_y, pix_x )[1] * 255.0;
	                    sv_pixel(2) = sv_image.at <cv::Vec3d> ( pix_y, pix_x )[0] * 255.0;

	                    /* push color vector */
	                    sv_color.push_back( sv_pixel );
                    }



                    /* Third picture match */
                    {
						pix_x += sv_bilinear_sample(sv_flow_23_u.pData, pix_x, pix_y, sv_width);
						pix_y += sv_bilinear_sample(sv_flow_23_v.pData, pix_x, pix_y, sv_width);
						pix_x = max(0.0, min((double)sv_width-1, pix_x));
						pix_y = max(0.0, min((double)sv_height-1, pix_y));
						sv_mat_3.push_back( sv_dense_geometry_cartesian( sv_width, sv_height, pix_x, pix_y ) );
                    }
                }
            }
        }
    }
#endif
/*
    source - scene methods
 */

    std::vector < Eigen::Vector3d > sv_dense_scene( std::vector < Eigen::Vector3d > const & sv_mat_1, std::vector < Eigen::Vector3d > const & sv_mat_2, std::vector < Eigen::Vector3d > const & sv_mat_3, Eigen::Vector3d const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3 ) {

        /* returned structure variable */
        std::vector < Eigen::Vector3d > sv_scene;

        /* parsing direction vectors */
        for ( long sv_parse( 0 ); sv_parse < sv_mat_1.size(); sv_parse ++ ) {

            /* compute and push optimised intersection */
            sv_scene.push_back( sv_dense_geometry_intersect( sv_mat_1[sv_parse], sv_mat_2[sv_parse], sv_mat_3[sv_parse], sv_cen_1, sv_cen_2, sv_cen_3 ) );

        }

        /* return computed scene */
        return( sv_scene );

    }

/*
    source - filtering methods
 */

    void sv_dense_filter( double const sv_tol, double const sv_max, std::vector < Eigen::Vector3d > const & sv_scene, std::vector < Eigen::Vector3i > const & sv_color, std::vector < Eigen::Vector3d > & sv_fscene, std::vector < Eigen::Vector3i > & sv_fcolor, std::vector < Eigen::Vector3d > const & sv_mat_1, std::vector < Eigen::Vector3d > const & sv_mat_2, std::vector < Eigen::Vector3d > const & sv_mat_3, Eigen::Vector3d const & sv_cen_1, Eigen::Vector3d const & sv_cen_2, Eigen::Vector3d const & sv_cen_3 ) {

        /* element radius variable */
        double sv_rad_1( 0.0 );
        double sv_rad_2( 0.0 );
        double sv_rad_3( 0.0 );

        /* element disparity variable */
        double sv_disp_1( 0.0 );
        double sv_disp_2( 0.0 );
        double sv_disp_3( 0.0 );

        /* parsing scene elements */
        for ( long sv_parse( 0 ); sv_parse < sv_scene.size(); sv_parse ++ ) {

            /* compute element radius */
            sv_rad_1 = sv_mat_1[sv_parse].transpose() * ( sv_scene[sv_parse] - sv_cen_1 );
            sv_rad_2 = sv_mat_2[sv_parse].transpose() * ( sv_scene[sv_parse] - sv_cen_2 );
            sv_rad_3 = sv_mat_3[sv_parse].transpose() * ( sv_scene[sv_parse] - sv_cen_3 );

            /* compute element disparity */
            sv_disp_1 = ( sv_cen_1 + ( sv_rad_1 * sv_mat_1[sv_parse] ) - sv_scene[sv_parse] ).norm();
            sv_disp_2 = ( sv_cen_2 + ( sv_rad_2 * sv_mat_2[sv_parse] ) - sv_scene[sv_parse] ).norm();
            sv_disp_3 = ( sv_cen_3 + ( sv_rad_3 * sv_mat_3[sv_parse] ) - sv_scene[sv_parse] ).norm();

            /* apply filtering condition */
            bool match = true;
            match &= sv_disp_1 <= sv_tol;
			match &= sv_disp_2 <= sv_tol;
			match &= sv_disp_3 <= sv_tol;
			match &= ( sv_rad_2 > 0.0 ) && ( sv_rad_2 < sv_max );

			if(sv_parse == 421447 || sv_parse == 420304){
				cout << match << endl;

			}

#ifdef BAD_MATCH_RED
			/* element selection - position */
			sv_fscene.push_back( sv_scene[sv_parse] );

			/* element selection - color */
			sv_fcolor.push_back( match ? sv_color[sv_parse] : Eigen::Vector3i(0,0,255));
#else
			if(match){
				/* element selection - position */
				sv_fscene.push_back( sv_scene[sv_parse] );

				/* element selection - color */
				sv_fcolor.push_back( sv_color[sv_parse] );
            }
#endif

        }

    }


#include <time.h>
void profile(string msg){
	static struct timespec specOld;
    struct timespec specNew;
    double t;

    clock_gettime(CLOCK_REALTIME, &specNew);

    if (specNew.tv_nsec >= 1000000000) {
    	specNew.tv_nsec -= 1000000000;
    	specNew.tv_sec++;
    }
    t = specNew.tv_sec - specOld.tv_sec + (specNew.tv_nsec - specOld.tv_nsec)*1e-9;
    std::cout << msg << " : " << t << endl;
    specOld = specNew;
}

/*
    source - main function
 */

//
//#include "NvOFCuda.h"
//#include "NvOFDataLoader.h"
//#include <memory>
//
//    int opticalFlow(){
//    	int gpuId = 0;
//        int nGpu = 0;
//        std::string inputFileName = "/home/dolu/pro/Optical_Flow_SDK_1.0.13/inputD";
//
//        CUDA_DRVAPI_CALL(cuInit(0));
//        CUDA_DRVAPI_CALL(cuDeviceGetCount(&nGpu));
//        if (gpuId < 0 || gpuId >= nGpu)
//        {
//            std::cout << "GPU ordinal out of range. Should be with in [" << 0 << ", " << nGpu - 1 << "]" << std::endl;
//            return 1;
//        }
//
//        CUdevice cuDevice = 0;
//        CUDA_DRVAPI_CALL(cuDeviceGet(&cuDevice, gpuId));
//        char szDeviceName[80];
//        CUDA_DRVAPI_CALL(cuDeviceGetName(szDeviceName, sizeof(szDeviceName), cuDevice));
//        std::cout << "GPU in use: " << szDeviceName << std::endl;
//
//        CUcontext cuContext = nullptr;
//        CUDA_DRVAPI_CALL(cuCtxCreate(&cuContext, 0, cuDevice));
//
//
//        std::unique_ptr<NvOFDataLoader> dataLoader = CreateDataloader(inputFileName);
//        uint32_t width = dataLoader->GetWidth();
//        uint32_t height = dataLoader->GetHeight();
//        uint32_t nFrameSize = width * height;
//
//        CUstream   inputStream = nullptr;
//        CUstream   outputStream = nullptr;
//
//        // Create Optical Flow object with desired frame width and height
//        NvOFObj nvOpticalFlow = NvOFCuda::Create(cuContext,
//           width,
//           height,
//           NV_OF_BUFFER_FORMAT_GRAYSCALE8,
//                                                    NV_OF_CUDA_BUFFER_TYPE_CUARRAY,
//                                                    NV_OF_CUDA_BUFFER_TYPE_CUARRAY,
//                                                    NV_OF_MODE_OPTICALFLOW,
//                                                    NV_OF_PERF_LEVEL_SLOW,
//                                                    inputStream,
//                                                    outputStream);
//
//
//        // Create input and output buffers
//        std::vector<NvOFBufferObj> inputBuffers = nvOpticalFlow->CreateBuffers(NV_OF_BUFFER_USAGE_INPUT, 2);
//        std::vector<NvOFBufferObj> outputBuffers = nvOpticalFlow->CreateBuffers(NV_OF_BUFFER_USAGE_OUTPUT, 1);
//
//
//        int curFrameIdx = 0;
//        inputBuffers[curFrameIdx]->UploadData(dataLoader->CurrentItem());
//        dataLoader->Next();
//        curFrameIdx ^= 1;
//
//        uint32_t nOutSize = outputBuffers[0]->getWidth() * outputBuffers[0]->getHeight();
//        std::unique_ptr<NV_OF_FLOW_VECTOR[]>pOut(new NV_OF_FLOW_VECTOR[nOutSize]);
//        if (pOut == nullptr)
//        {
//            std::ostringstream err;
//            err << "Failed to allocate output host memory of size " << nOutSize * sizeof(NV_OF_FLOW_VECTOR) << " bytes" << std::endl;
//            throw std::bad_alloc();
//        }
//
//        // Read frames and calculate optical flow vectors on consecutive frames
//        for (; !dataLoader->IsDone(); dataLoader->Next())
//        {
//            inputBuffers[curFrameIdx]->UploadData(dataLoader->CurrentItem());
//
//            nvOpticalFlow->Execute(inputBuffers[curFrameIdx ^ 1].get(),
//                                                inputBuffers[curFrameIdx].get(),
//                                                outputBuffers[0].get());
//
//            outputBuffers[0]->DownloadData(pOut.get());
//
//            // Use the flow vectors generated in buffer pOut.get() for further
//            // processing or save to a file.
//
//            curFrameIdx = curFrameIdx ^ 1;
//        }
//
//
//        CUDA_DRVAPI_CALL(cuCtxDestroy(cuContext));
//
//        std::cout << "MIAOU" << endl;
//        return 0;
//    }


    int main( int argc, char ** argv ) {
		profile("init");
    	//exit(opticalFlow());


        /* matrix variable */
        cv::Mat sv_img_prev;
        cv::Mat sv_img_middle;
        cv::Mat sv_img_next;

        /* matrix variable */
        cv::Mat sv_mask;

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
        std::vector < Eigen::Vector3d > sv_fscene;

        /* color variable */
        std::vector < Eigen::Vector3i > sv_color;
        std::vector < Eigen::Vector3i > sv_fcolor;

        /* filtering variable */
        double sv_tol( 0.0 );
        double sv_max( 0.0 );


        /* check consistency */
#ifdef OPTICAL_FLOW_CPP
        if ( argc != 8 ) {

            /* display message */
            std::cerr << "scanvan : error : wrong usage" << std::endl;

            /* display quick help */
            std::cerr << "./Dense [image 1 path] [image 2 path] [image 3 path] [pose estimation file] [scene export path] [mask image path] [image scale]" << std::endl;

            /* send message */
            return( 1 );

        }
#endif

#ifdef OPTICAL_FLOW_NVIDIA
        if ( argc != 10 ) {

            /* display message */
            std::cerr << "scanvan : error : wrong usage" << std::endl;

            /* display quick help */
            std::cerr << "./Dense [image 1 path] [image 2 path] [image 3 path] [pose estimation file] [scene export path] [mask image path] [image scale] [nvidia flow 12] [nvidia flow 23]" << std::endl;

            std::cerr << "got " << argc << " arguments = "   << std::endl << argv[1] << std::endl;
            /* send message */
            return( 1 );

        }
#endif
        float sv_scale = atof( argv[7] );

        /* import estimation parameters */
        sv_dense_io_pose( argv[4], sv_r12, sv_t12, sv_r23, sv_t23 );

        /* import image */
        sv_img_prev   = sv_dense_io_image( argv[1], sv_scale );
        sv_img_middle = sv_dense_io_image( argv[2], sv_scale );
        sv_img_next   = sv_dense_io_image( argv[3], sv_scale );

        /* extract image reference */
        sv_img_width  = sv_img_middle.cols;
        sv_img_height = sv_img_middle.rows;
        sv_img_depth  = sv_img_middle.channels();

        /* check image */
        sv_dense_consistent_image( sv_img_prev, sv_img_width, sv_img_height, sv_img_depth );
        sv_dense_consistent_image( sv_img_next, sv_img_width, sv_img_height, sv_img_depth );

        /* import mask */
        sv_mask = sv_dense_io_mask( argv[6], sv_scale );


        profile("args");
    # pragma omp parallel sections
    {

    # pragma omp section
    {
        /* compute optical flow : image 2 -> 1 */
#ifdef OPTICAL_FLOW_CPP
        sv_dense_flow( sv_img_middle, sv_img_prev, sv_img_width, sv_img_height, sv_img_depth, sv_flow_21_u, sv_flow_21_v );
#endif
#ifdef OPTICAL_FLOW_NVIDIA
        sv_dense_flow_nvidia(argv[8], sv_scale, sv_flow_21_u, sv_flow_21_v);
#endif
    }

    # pragma omp section
    {

        /* compute optical flow : image 2 -> 3 */
#ifdef OPTICAL_FLOW_CPP
        sv_dense_flow( sv_img_middle, sv_img_next, sv_img_width, sv_img_height, sv_img_depth, sv_flow_23_u, sv_flow_23_v );
#endif
#ifdef OPTICAL_FLOW_NVIDIA
        sv_dense_flow_nvidia(argv[9], sv_scale, sv_flow_23_u, sv_flow_23_v);
#endif
    }

    } /* parallel sections */
    	profile("sv_dense_flow_nvidia");

        /* compute matches */
        sv_dense_match( sv_img_middle, sv_mask, sv_img_width, sv_img_height, sv_flow_21_u, sv_flow_21_v, sv_flow_23_u, sv_flow_23_v, sv_mat_1, sv_mat_2, sv_mat_3, sv_color );
        profile("sv_dense_match");

        /* compute common frame - aligned on first camera */
        sv_dense_geometry_common( sv_mat_1, sv_mat_2, sv_mat_3, sv_cen_1, sv_cen_2, sv_cen_3, sv_r12, sv_t12, sv_r23, sv_t23 );
        profile("sv_dense_geometry_common");

        /* compute scene */
        sv_scene = sv_dense_scene( sv_mat_1, sv_mat_2, sv_mat_3, sv_cen_1, sv_cen_2, sv_cen_3 );
        profile("sv_dense_scene");

        /* compute filtering tolerence values */
        //sv_tol = sv_t12.norm() + sv_t23.norm();
        sv_tol = sv_dense_geometry_amplitude( sv_cen_1, sv_cen_2, sv_cen_3 );
        sv_max = sv_tol *  20.0;
        sv_tol = sv_tol / 150.0;
        profile("sv_dense_geometry_amplitude");

        /* filter scene */
        sv_dense_filter( sv_tol, sv_max, sv_scene, sv_color, sv_fscene, sv_fcolor, sv_mat_1, sv_mat_2, sv_mat_3, sv_cen_1, sv_cen_2, sv_cen_3 );
        profile("sv_dense_filter");

        /* export computed scene */
        sv_dense_io_scene( argv[5], sv_fscene, sv_fcolor );
        profile("sv_dense_io_scene");

        /* send message */
        return( 0 );

    }

