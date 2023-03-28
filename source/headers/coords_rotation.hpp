//  This header file contains functions and classes for working with spherical and Cartesian coordinates/
//
//  1) ALPHA_DEF, BETA_DEF, GAMMA_DEF -- Euler angles in degrees
//  2) eps                            -- offset value
//  3) Euler_rotation_info            -- class that stores information about rotation and performs it
//  4) Spherical_Coords               -- class that stores spherical coords
//  5) Cartesian_Coords               -- class that stores Cartesian coords
//  6) Overloaded operators <<        -- functions that prints spherical and Cartesian coordinates to output stream
//  7) Overloaded operator ==         -- checks equality of spherical coordinates
//  8) From_spherical_to_Cartesian    -- function that performs transformation from spherical to Cartesian coordinates (with unit radius)
//  9) From_Cartesian_to_spherical    -- function that performs transformation from Cartesian to spherical coordinates (with unit radius)
// 10) Rotate_Cartesian               -- function that performs rotation of point with Cartesian coordinates
// 11) Rotate_Spherical               -- function that performs rotation of point with spherical coordinates

#pragma once
#include <vector>
#include <cmath>
#include <ostream>
#include <proj_api.h>
#include "inmost.h"
#include "defines.hpp"

const double ALPHA_DEF               = -30.0;  // default alpha for Euler rotation in degrees
const double BETA_DEF                = -90.0;  // default beta for Euler rotation in degrees
const double GAMMA_DEF               = 0.0;    // default gamma for Euler rotation in degrees
const double eps                     = 1e-10;  // offset to avoid zero values
const int    OUTPUT_DOUBLE_PRECISION = 13;     // total number of digits in output double

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  #
// #    !!! Purpose of class - store the information about rotation !!!                                  #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  #
// #                                                                                                     #
// #    Constructor input : Euler angles: alpha, beta and gamma                                          #
// #    Contains:                                                                                        #
// #    ALPHA   - real rotation angle about z-axis in degrees can be obtained by call .Get_ALPHA()       #
// #    BETA    - real rotation angle about new y-axis in degrees can be obtained by call .Get_BETA()    #
// #    GAMMA   - real rotation angle about new z-axis in degrees can be obtained by call .Get_GAMMA()   #
// #    FORWARD - 3x3 real matrix that perform rotation from model to geographical Cartesian coordinates #
// #          (appears after a call constructor) can be obtained by call .Get_FORWARD()                  #
// #    REVERSE - 3x3 real matrix that perform rotation from geographical to model Cartesian coordinates #
// #          (appears after a call constructor) can be obtained by call .Get_REVERSE()                  #
// #                                                                                                     #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
namespace SIMUG
{

template<typename RealType>
class Euler_rotation_info
{
public:
  Euler_rotation_info(RealType alpha, RealType beta, RealType gamma)
  : ALPHA(alpha)
  , BETA (beta)
  , GAMMA(gamma)
  {
    this->GenerateRotationMatricies();
  }

  const std::vector<std::vector<RealType>>& Get_FORWARD() const
   {
     return FORWARD;
   }

  const std::vector<std::vector<RealType>>& Get_REVERSE() const
   {
     return REVERSE;
   }

   RealType Get_ALPHA() const
   {
     return ALPHA;
   }

   RealType Get_BETA() const
   {
     return BETA;
   }

   RealType Get_GAMMA() const
   {
     return BETA;
   }
private:
 const RealType ALPHA;                        // rotation about z-axis
 const RealType BETA;                         // rotation about new y-axis
 const RealType GAMMA;                        // rotation about new z-axis
 std::vector<std::vector<RealType>> FORWARD;  // rotation matricies to geographic coordinates
 std::vector<std::vector<RealType>> REVERSE;  // rotation matricies to model coordinates

 void GenerateRotationMatricies()
 {
   RealType ALPHA_RAD = ALPHA*M_PI/180.0;
   RealType BETA_RAD = BETA*M_PI/180.0;
   RealType GAMMA_RAD = GAMMA*M_PI/180.0;

   FORWARD.push_back({
     cos(GAMMA_RAD)*cos(BETA_RAD)*cos(ALPHA_RAD) - sin(GAMMA_RAD)*sin(ALPHA_RAD),
     -sin(GAMMA_RAD)*cos(BETA_RAD)*cos(ALPHA_RAD) - cos(GAMMA_RAD)*sin(ALPHA_RAD),
     sin(BETA_RAD)*cos(ALPHA_RAD)
   });

   FORWARD.push_back({
     cos(GAMMA_RAD)*cos(BETA_RAD)*sin(ALPHA_RAD) + sin(GAMMA_RAD)*cos(ALPHA_RAD),
     -sin(GAMMA_RAD)*cos(BETA_RAD)*sin(ALPHA_RAD) + cos(GAMMA_RAD)*cos(ALPHA_RAD),
     sin(BETA_RAD)*sin(ALPHA_RAD)
   });

   FORWARD.push_back({
     -cos(GAMMA_RAD)*sin(BETA_RAD),
     sin(GAMMA_RAD)*sin(BETA_RAD),
     cos(BETA_RAD)
   });

   REVERSE.push_back({
     cos(ALPHA_RAD)*cos(BETA_RAD)*cos(GAMMA_RAD) - sin(ALPHA_RAD)*sin(GAMMA_RAD),
     sin(ALPHA_RAD)*cos(BETA_RAD)*cos(GAMMA_RAD) + cos(ALPHA_RAD)*sin(GAMMA_RAD),
     -sin(BETA_RAD)*cos(GAMMA_RAD)
   });

   REVERSE.push_back({
     -cos(ALPHA_RAD)*cos(BETA_RAD)*sin(GAMMA_RAD) - sin(ALPHA_RAD)*cos(GAMMA_RAD),
     -sin(ALPHA_RAD)*cos(BETA_RAD)*sin(GAMMA_RAD) + cos(ALPHA_RAD)*cos(GAMMA_RAD),
     sin(BETA_RAD)*sin(GAMMA_RAD)
   });

   REVERSE.push_back({
     cos(ALPHA_RAD)*sin(BETA_RAD),
     sin(ALPHA_RAD)*sin(BETA_RAD),
     cos(BETA_RAD)
   });

   FORWARD.shrink_to_fit();
   REVERSE.shrink_to_fit();
 }
};


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                            #
// #    !!! Purpose of class - store spherical coordinates of point !!!                            #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                            #
// #                                                                                               #
// #     Contains:                                                                                 #
// #     x_      - real first coordinate  (longitude) in degrees can be obtained by call .Get_x()  #
// #     y_      - real second coordinate (latitude) in degrees can be obtained by call .Get_y()   #
// #     +=      - increment real number to both components                                        #
// #     -=      - decrement real number to both components                                        #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
template<typename RealType>
class Spherical_Coords
{
public:
  Spherical_Coords()
  {}

  Spherical_Coords(RealType x, RealType y)
  : x_(x)
  , y_(y)
  {}

  RealType Get_x() const
  {
    return x_;
  }

  RealType Get_y() const
  {
    return y_;
  }

  void operator += (RealType increment)
  {
    x_ += increment;
    y_ += increment;
  }

  void operator -= (RealType decrement)
  {
    x_ -= decrement;
    y_ -= decrement;
  }

  void SetWater()
  {
    is_water = true;
  }

  bool IsWater() const
  {
    return is_water;
  }

private:
  RealType x_;
  RealType y_;
  bool is_water = false;
};



// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #    !!! Purpose of class - store Cartesian coordinates of point !!!      #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                         #
// #     Contains:                                                           #
// #     x_      - real first coordinate can be obtained by call .Get_x()    #
// #     y_      - real second coordinate can be obtained by call .Get_y()   #
// #     z_      - real third coordinate can be obtained by call .Get_z()    #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
template<typename RealType>
class Cartesian_Coords
{
public:
  Cartesian_Coords()
  {}

  Cartesian_Coords(RealType x, RealType y, RealType z)
  : x_(x)
  , y_(y)
  , z_(z)
  {}

  RealType Get_x() const
  {
    return x_;
  }

  RealType Get_y() const
  {
    return y_;
  }

  RealType Get_z() const
  {
    return z_;
  }

private:
  RealType x_;
  RealType y_;
  RealType z_;
};

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - print spherical coords to output stream (cout or file) !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input : point with spherical coordinates                                                   #
// #    output: reference to output stream                                                         #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
template<typename RealType>
std::ostream& operator<<(std::ostream& stream, const Spherical_Coords<RealType>& cords)
{
  stream.precision(OUTPUT_DOUBLE_PRECISION);
  stream << cords.Get_x() << ' ' << cords.Get_y(); //<< cords.IsWater() << ' ';
  return stream;
}

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - print Cartesian coords to output stream (cout or file) !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input : point with Cartesian coordinates                                                   #
// #    output: reference to output stream                                                         #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
template<typename RealType>
std::ostream& operator<<(std::ostream& stream, const Cartesian_Coords<RealType>& cords)
{
  stream.precision(OUTPUT_DOUBLE_PRECISION);
  stream << cords.Get_x() << ' ' << cords.Get_y() << ' ' << cords.Get_z();
  return stream;
}

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - identify equality of two spherical points              !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input1 : lhs - first spherical point                                                       #
// #    input2 : rhs - second spherical point                                                      #
// #    output: boolean (true for =, false for !=)                                                 #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
template<typename RealType>
bool operator==(const Spherical_Coords<RealType>& lhs, const Spherical_Coords<RealType>& rhs)
{
  return ((lhs.Get_x() == rhs.Get_x()) and (lhs.Get_y() == rhs.Get_y()));
}

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - transform spherical point to Cartesian one             !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input : Sph_coords - point with spherical coords                                           #
// #    output: point with Cartesian coords                                                        #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
template<typename RealType>
Cartesian_Coords<RealType> From_spherical_to_Cartesian(const Spherical_Coords<RealType>& Sph_coords)
{
  double x_rad = Sph_coords.Get_x()*M_PI/180.0;
  double y_rad = Sph_coords.Get_y()*M_PI/180.0;

  return {cos(y_rad)*cos(x_rad),
          cos(y_rad)*sin(x_rad),
          sin(y_rad)};
}


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - transform Cartesian point to spherical one             !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input : Crt_coords - point with Cartesian coords                                           #
// #    output: point with spherical coords                                                        #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
template<typename RealType>
Spherical_Coords<RealType> From_Cartesian_to_spherical(const Cartesian_Coords<RealType>& Crt_coords)
{
  RealType xx = Crt_coords.Get_x();
  RealType yy = Crt_coords.Get_y();
  RealType zz = Crt_coords.Get_z();

  RealType phi, theta, costn;

  theta = asin(zz)*180/M_PI;
  costn = sqrt(1.0 - zz*zz);

  if ((xx > 0.0) and (yy > 0.0))
    {
      if (xx < yy)
      {
        phi = acos(xx/costn)*180.0/M_PI;
      }
      else
      {
        phi = asin(yy/costn)*180.0/M_PI;
      }
    }
    else if ((xx < 0.0) and (yy > 0.0))
    {
      if (fabs(xx) < yy)
      {
        phi = 180.0 - acos(fabs(xx)/costn)*180.0/M_PI;
      }
      else
      {
        phi = 180.0 - asin(yy/costn)*180.0/M_PI;
      }
    }
    else if ((xx < 0.0) and (yy < 0.0))
    {
      if (fabs(xx) < fabs(yy))
      {
        phi = -180.0 + acos(fabs(xx)/costn)*180/M_PI;
      }
      else
      {
        phi = -180.0 + asin(fabs(yy)/costn)*180/M_PI;
      }
    }
    else if ((xx > 0.0) and (yy < 0.0))
    {
      if (xx < fabs(yy))
      {
        phi = -acos(fabs(xx)/costn)*180.0/M_PI;
      }
      else
      {
        phi = -asin(fabs(yy)/costn)*180.0/M_PI;
      }
    }

    //if (phi < 0.0)
    //{
    //  phi += 360.0;
    //}

    return {phi, theta};
}


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - calculate new Cartesian coordinates of point after     !!!      #
// #     !!!                       rotation with given rotation matrix                    !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input1: cartesian_old   - Cartesian coordinates of point before rotation                   #
// #    input2: rot             - 3x3 rotation matrix                                              #
// #    output: Cartesian coordinate of point after rotation                                       #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
template<typename RealType>
Cartesian_Coords<RealType> Rotate_Cartesian(const Cartesian_Coords<RealType>&  cartesian_old, const std::vector<std::vector<RealType>>& rot)
{
  return {rot[0][0]*cartesian_old.Get_x() + rot[0][1]*cartesian_old.Get_y() + rot[0][2]*cartesian_old.Get_z(),
          rot[1][0]*cartesian_old.Get_x() + rot[1][1]*cartesian_old.Get_y() + rot[1][2]*cartesian_old.Get_z(),
          rot[2][0]*cartesian_old.Get_x() + rot[2][1]*cartesian_old.Get_y() + rot[2][2]*cartesian_old.Get_z()
  };
}

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - calculate new spherical coordinates of point after     !!!      #
// #     !!!                       rotation with given rotation matrix                    !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input1: spherical_old   - spherical coordinates of point before rotation                   #
// #    input2: rot             - 3x3 rotation matrix                                              #
// #    output: spherical coordinate of point after rotation                                       #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
template<typename RealType>
Spherical_Coords<RealType> Rotate_Spherical(const Spherical_Coords<RealType>& spherical_old, const std::vector<std::vector<RealType>>& rot)
{

  Spherical_Coords<RealType>  spherical_new;   //output
  Spherical_Coords<RealType>  spherical_tmp;   //temporary
  Cartesian_Coords<RealType>  cartesian_old;   //old Cartesian coordinate
  Cartesian_Coords<RealType>  cartesian_new;   //new Cartesian coordinate

  spherical_tmp = spherical_old;

  //avoid trouble of an exactly zero angle by adding offset
  spherical_tmp += eps;

  //spherical to Cartesian
  cartesian_old = From_spherical_to_Cartesian(spherical_tmp);

  //new Cartesian coordinates after rotation
  cartesian_new = Rotate_Cartesian(cartesian_old, rot);

  //Cartesian to spherical
  spherical_new = From_Cartesian_to_spherical(cartesian_new);

   //avoid trouble of an exactly zero angle by subtracting offset
  spherical_new -= eps;

  return spherical_new;
}

template<typename RealType>
std::vector<RealType> from_model_2_geo(RealType model_lon, RealType model_lat)
{
  Euler_rotation_info<RealType> rotation(ALPHA_DEF, BETA_DEF, GAMMA_DEF);
  RealType lat, lon;
  Spherical_Coords<RealType> model_forw(model_lon, model_lat);
  Spherical_Coords<RealType> geo_forw = Rotate_Spherical<RealType>(model_forw, rotation.Get_FORWARD());

  lon = geo_forw.Get_x();
  lat = geo_forw.Get_y();

  return {lon, lat};
}


template<typename RealType>
std::vector<RealType> from_geo_2_model(RealType lon, RealType lat)
{
  Euler_rotation_info<RealType> rotation(ALPHA_DEF, BETA_DEF, GAMMA_DEF);
  RealType model_lat, model_lon;
  Spherical_Coords<RealType> geo_rev(lon, lat);
  Spherical_Coords<RealType> model_rev = Rotate_Spherical<RealType>(geo_rev, rotation.Get_REVERSE());

  model_lon = model_rev.Get_x();
  model_lat = model_rev.Get_y();

  return {model_lon, model_lat};
}


template<typename RealType>
std::vector<RealType> from_geo_2_model_vec(RealType geo_vec_lon, RealType geo_vec_lat,
                                           RealType geo_lon, RealType geo_lat)
{
  std::vector<RealType> m = from_geo_2_model(geo_lon, geo_lat);
  RealType LON_MOD = m[0];
  RealType LAT_MOD = m[1];

  RealType X_GEO, Y_GEO, Z_GEO;
  RealType radian = M_PI/180.0;
  Euler_rotation_info<RealType> rotation(ALPHA_DEF, BETA_DEF, GAMMA_DEF);

  //convert geogr. long/lat into geogr. cartesian
  X_GEO = -geo_vec_lon*sin(geo_lon*radian) -
  geo_vec_lat*cos(geo_lon*radian)*sin(geo_lat*radian);

  Y_GEO = geo_vec_lon*cos(geo_lon*radian) -
  geo_vec_lat*sin(geo_lon*radian)*sin(geo_lat*radian);
 
  Z_GEO = geo_vec_lat*cos(geo_lat*radian);

  //rotate geogr. cartesian into model cartesian
  Cartesian_Coords<RealType> c(X_GEO, Y_GEO, Z_GEO);
  Cartesian_Coords<RealType> rc = Rotate_Cartesian<RealType>(c, rotation.Get_REVERSE());
  
  
  double FIELD_X, FIELD_Y;
  double X_MOD = rc.Get_x();
  double Y_MOD = rc.Get_y();
  double Z_MOD = rc.Get_z();

  //convert model cartesian into model long/lat
  FIELD_X = - X_MOD*sin(LON_MOD*radian)
            + Y_MOD*cos(LON_MOD*radian);
            
  FIELD_Y = - X_MOD*cos(LON_MOD*radian)*sin(LAT_MOD*radian)
            - Y_MOD*sin(LON_MOD*radian)*sin(LAT_MOD*radian)
            + Z_MOD*cos(LAT_MOD*radian);

  return {FIELD_X, FIELD_Y};
}

std::vector<double> from_geo_2_topaz(double lon, double lat);
std::vector<double> from_topaz_2_geo(double x, double y);
}