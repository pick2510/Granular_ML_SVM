#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

enum calc_layout
{
    rad1,
    rad2,
    Radstar,
    Penor,
    Petan,
    Ftan,
    Fnor,
    kn,
    kt,
    Petot,
};

enum ContactTXTColumns
{
    p1_x,
    p1_y,
    p1_z,
    p2_x,
    p2_y,
    p2_z,
    p1_vx, 
    p1_vy,
    p1_vz,
    p2_vx,
    p2_vy,
    p2_vz,
    p1_id,
    p2_id,
    is_periodic,
    c_force_x,
    c_force_y,
    c_force_z, 
    cn_force_x,
    cn_force_y,
    cn_force_z,
    ct_force_x,
    ct_force_y,
    ct_force_z,
    c_torque_x,
    c_torque_y,
    c_torque_z,
    disp_x,
    disp_y,
    disp_z,
    contact_area,
    contact_overlap,
    sliding_contact, 
};

template <class Container>
void split(const std::string &str, Container &cont)
{
    std::istringstream iss(str);
    std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(),
              std::back_inserter(cont));
}

#endif