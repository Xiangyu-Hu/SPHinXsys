
#include "utilities.h"

namespace SPH
{
    Vecd rotatePointAroundPoint(Vecd point, Real angle, Vecd rotate_center)
    {
        Real cos_value = std::cos(angle);
        Real sin_value = std::sin(angle);
        Vecd temp = point - rotate_center;
        Real x_new = temp[0] * cos_value - temp[1] * sin_value;
        Real y_new = temp[0] * sin_value + temp[1] * cos_value;

        temp = Vecd(x_new + rotate_center[0], y_new + rotate_center[1], 0);
        return temp;
    }

    BoundingBox calculateNewBoundingBox(BoundingBox old_box, Real rotate_angle, Vecd rotate_center)
    {
        Vecd &old_lower_left = old_box.first;
        Vecd &oldTopRight = old_box.second;
        Real delta_x = oldTopRight[0] - old_lower_left[0];
        Real delta_y = oldTopRight[1] - old_lower_left[1];
        Vecd old_top_left = old_lower_left + Vecd(0, delta_y);
        Vecd old_lower_right = old_lower_left + Vecd(delta_x, 0);

        auto p1 = rotatePointAroundPoint(old_lower_left, rotate_angle, rotate_center);
        auto p2 = rotatePointAroundPoint(oldTopRight, rotate_angle, rotate_center);
        auto p3 = rotatePointAroundPoint(old_top_left, rotate_angle, rotate_center);
        auto p4 = rotatePointAroundPoint(old_lower_right, rotate_angle, rotate_center);

        Vecd new_lower_left(SMIN(p1[0], p2[0], p3[0], p4[0]), SMIN(p1[1], p2[1], p3[1], p4[1]));
        Vecd new_top_right(SMAX(p1[0], p2[0], p3[0], p4[0]), SMAX(p1[1], p2[1], p3[1], p4[1]));

        return std::make_pair(new_lower_left, new_top_right);
    }
}
