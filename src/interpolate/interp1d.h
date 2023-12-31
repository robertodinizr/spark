#pragma once

#include "containers/containers.h"

#include <cstring>
#include <utility>
#include <algorithm>
#include <limits>

namespace sci::interpolate
{

    enum class interp1d_type
    {
        Nearest,
        NearestUp,
        Linear,
        Next
    };

    template <typename T>
    class interp1d
    {

    public:
        interp1d(sci::buffer<T> &&x, sci::buffer<T> &&y, interp1d_type type) : xf(std::move(x)), yf(std::move(y)), typef(type) {}

        T operator()(T xp)
        {
            switch (typef)
            {
            case interp1d_type::Nearest:
                return nearest(xp);
                break;
            
            case interp1d_type::NearestUp:
                return nearest_up(xp);
                break;

            default:
                break;
            }

            return std::numeric_limits<T>::quiet_NaN();
        }

    private:
        T nearest(T xp)
        {
            if (xp >= xf.data[yf.size - 1])
                return yf.data[yf.size - 1];

            if (xp < xf.data[0])
                return std::numeric_limits<T>::quiet_NaN();

            auto lower = std::upper_bound(xf.begin(), xf.end(), xp);
            auto idx = std::distance(xf.begin(), lower) - 1;
            return yf.data[idx];
        }

        T nearest_up(T xp)
        {
            if (xp > xf.data[yf.size - 1])
                return std::numeric_limits<T>::quiet_NaN();

            if (xp <= xf.data[0])
                return xf.data[0];

            auto lower = std::lower_bound(xf.begin(), xf.end(), xp);
            auto idx = std::distance(xf.begin(), lower);
            return yf.data[idx];
        }

        sci::buffer<T> xf, yf;
        interp1d_type typef;
    };
}