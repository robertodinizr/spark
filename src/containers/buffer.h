#pragma once

#include <cstddef>
#include <cstring>
#include <initializer_list>

namespace sci::containers
{
    template <typename T>
    struct buffer
    {
        T *data = nullptr;
        size_t size = 0;

        explicit buffer() {}
        buffer(T* _data, size_t _size) : buffer(_size)
        {
            for(size_t i = 0; i < size; i++)
                data[i] = _data[i];
        }
        
        buffer(const std::initializer_list<T>& vec) : buffer((T*) vec.begin(), vec.size()){}
        explicit buffer(size_t _size) : data(new T[_size]), size(_size) {}


        ~buffer()
        {
            delete data;
        }

        // Copy constructor and assignment operator deleted
        buffer(const buffer &) = delete;
        buffer &operator=(const buffer &) = delete;

        // Move constructor and assignment operator
        buffer(buffer &&other) noexcept : data(other.data), size(other.size)
        {
            other.data = nullptr;
            other.size = 0;
        }

        buffer &operator=(buffer &&other) noexcept
        {
            if (this != &other)
            {
                delete data;
                data = other.data;
                size = other.size;
                other.data = nullptr;
                other.size = 0;
            }
            return *this;
        }

        buffer copy()
        {
            return buffer(data, size);
        }

        T *begin() const { return data; }
        T *end() const { return data + size; }
    };
}