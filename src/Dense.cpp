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
        if ( argc != 3 ) {

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

        /* import and check image */
        if ( ! ( sv_img_middle = cv::imread( argv[2], CV_LOAD_IMAGE_COLOR ) ).data ) {

            /* display message */
            std::cerr << "scanvan : error : unable to import image" << std::endl;

            /* send message */
            return( 1 );

        }

        /* import and check image */
        if ( ! ( sv_img_next = cv::imread( argv[3], CV_LOAD_IMAGE_COLOR ) ).data ) {

            /* display message */
            std::cerr << "scanvan : error : unable to import image" << std::endl;

            /* send message */
            return( 1 );

        }

        /* send message */
        return( 0 );

    }
