/** 
 *  \file
 *  Color features
 * ----------------------------------------------------------------------------
 *  \date       2016
 *  \author     Simon Rohou
 *  \copyright  Copyright 2019 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef __TUBEX_COLORS_H__
#define __TUBEX_COLORS_H__

#include <string>

namespace tubex
{
  /**
   * \struct rgb
   * \brief Represents an RGB value
   */
  typedef struct
  {
    float r;      ///< red, value between 0 and 1
    float g;      ///< green, value between 0 and 1
    float b;      ///< blue, value between 0 and 1
    float alpha;  ///< opacity, value between 0 (transparent) and 1 (opaque)
  } rgb;

  /**
   * \struct hsv
   * \brief Represents an HSV value
   */
  typedef struct
  {
    float h;      ///< hue, angle value in degrees between 0 and 359
    float s;      ///< saturation, a value between 0 and 1
    float v;      ///< value (lightness), a value between 0 and 1
    float alpha;  ///< opacity, value between 0 (transparent) and 1 (opaque)
  } hsv;

  /**
   * \brief Represents an RGB value in a HTML standard
   *
   * \param rgb_value
   * \param prefix optional caracters ("#" by default)
   * \return the HTML string
   */
  const std::string rgb2hex(rgb rgb_value, const char* prefix = "#");

  /**
   * \brief Converts HSV to RGB
   *
   * \param hsv_value
   * \return RGB value
   */
  rgb hsv2rgb(hsv hsv_value);

  /**
   * \brief Converts RGB to HSV
   *
   * \param rgb_value
   * \return HSV value
   */
  hsv rgb2hsv(rgb rgb_value);

  /**
   * \brief Makes an RGV value from integers
   *
   * \param r red value, integer between 0 and 255
   * \param g green value, integer between 0 and 255
   * \param b blue value, integer between 0 and 255
   * \param alpha opacity value, integer between 0 (transparent) and 255 (opaque)
   * \return RGB value
   */
  rgb make_rgb(int r, int g, int b, int alpha = 255);

  /**
   * \brief Makes an RGV value from floats
   *
   * \param r red value, float between 0. and 1.
   * \param g green value, float between 0. and 1.
   * \param b blue value, float between 0. and 1.
   * \param alpha opacity value, float between 0. (transparent) and 1. (opaque)
   * \return RGB value
   */
  rgb make_rgb(float r, float g, float b, float alpha = 1.);

  /**
   * \brief Returns the RGB value of the Haxby color map corresponding to some ratio
   *
   * \param ratio a value between 0 and 1
   * \return RGB value in Haxby color map
   */
  rgb haxby_color_map(float ratio);

  /**
   * \brief Returns the RGB value of the DEM color map corresponding to some ratio
   *
   * \param ratio a value between 0 and 1
   * \return RGB value in DEM color map
   */
  rgb dem_color_map(float ratio);
}

#endif