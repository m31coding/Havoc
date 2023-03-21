#ifndef MY_ARRAY
#define MY_ARRAY

#include "assert.h"
#include <cstdio>

template<class T>
class myArray
{
private:

    T* m_data; ///< the array to store
    unsigned int m_size; ///< the size of the array
    unsigned int m_counter; ///< temporary position, array content with an index smaller then counter is possibly in use
    unsigned int m_deactivated; ///< number of elements with an index smaller then counter and not in use

public:

    myArray(unsigned int array_size)
    {
        m_size = array_size;
        m_counter = 0;
        m_deactivated = 0;
        m_data = new T[array_size];
    }

    myArray()
    {}

    ~myArray()
    {

        delete[] m_data;
    }

    /// delete a single element
    inline void del(T* element)
    {
        assert(m_counter - 1 > 0 && m_counter - 1 < m_size);
        *element = m_data[m_counter - 1];
        m_counter--;
    }

    /// increase the counter
    inline void counterInc()
    {
        m_counter++;
    }

    /// increase the deactivated counter
    inline void deactivatedInc()
    {
        m_deactivated++;
    }

    /// reset the counter and the number of deactivated elements
    inline void reset()
    {
        m_counter = 0;
        m_deactivated = 0;
    }

    /// get the number of active elements
    inline unsigned int NOFactive()
    {
        return (m_counter - m_deactivated);
    }

    /// get the counter
    inline unsigned int counter()
    {
        return m_counter;
    }

    /// get the maximum size of the array
    inline unsigned int max_size()
    {
        return m_size;
    }

    /// get the element with index counter
    inline T& current()
    {
        assert(m_counter < m_size); //check range

        return m_data[m_counter];
    }

    /// get the element with the given index
    inline T& at(unsigned int i)
    {
        assert(i < m_size);
        return m_data[i];
    }

    /// operator[]
    inline operator const T*(void) const
    {
        return m_data;
    }

    inline operator T*(void)
    {
        return m_data;
    }

    /// get number of available elements
    inline int free()
    {
        return (m_size - m_counter);
    }
};

#endif