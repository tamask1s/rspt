/***************************************************************************
* Copyright 2020 Tamas Levente Kis                                         *
*                                                                          *
* Licensed under the Apache License, Version 2.0 (the "License");          *
* you may not use this file except in compliance with the License.         *
* You may obtain a copy of the License at                                  *
*                                                                          *
*     http://www.apache.org/licenses/LICENSE-2.0                           *
*                                                                          *
* Unless required by applicable law or agreed to in writing, software      *
* distributed under the License is distributed on an "AS IS" BASIS,        *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
* See the License for the specific language governing permissions and      *
* limitations under the License.                                           *
***************************************************************************/

#ifndef _ZAX_TENSOR_H_
#define _ZAX_TENSOR_H_

/** containers storing their data in continuous memory fields. */

#define OPERATOR_EQ()\
public:\
    virtual void operator = (const char* a_json) {\
        zax_from_json(a_json);\
    }\
    virtual void operator = (const std::string& a_json) {\
        zax_from_json(a_json.c_str());\
    }\

template<typename T>
class tensor_base
{
public:
    int d1 = 0;
    int d2 = 0;
    int d3 = 0;
    int d4 = 0;

    T* d1d = 0;
    T** d2d = 0;
    T*** d3d = 0;
    T**** d4d = 0;

    virtual T* data() const = 0;

    virtual void detach_data() = 0;

    virtual std::vector<int> shape() const = 0;

    virtual bool reshape(const std::vector<int>& a_sizes) = 0;

    virtual void resize(const std::vector<int>& a_sizes, T* a_data = 0, int a_wrap_around_bytes = false, bool a_reshaping = false) = 0;

    virtual void zax_from_json(const char* a_json, std::string* a_err_stream = 0) = 0;

    virtual int zax_to_json(char* a_json, int a_alloc_size, int a_deep) const = 0;

    virtual void fill(const T& a_val) = 0;

    tensor_base()
    {}

    virtual ~tensor_base()
    {}

    virtual std::string shape_s() const
    {
        std::vector<int> sh = shape();
        std::string res = "{";
        char tmp[16];
        if (!sh.empty())
        {
            for (unsigned int i = 0; i < sh.size() - 1; ++i)
            {
                sprintf(tmp, "%d,", sh[i]);
                res += tmp;
            }
            sprintf(tmp, "%d", sh.back());
            res += tmp;
        }
        res += "}";
        return res;
    }

    virtual int numel() const
    {
        std::vector<int> s = shape();
        return std::accumulate(s.begin(), s.end(), 1, std::multiplies<int>());
    }

    virtual bool operator == (const tensor_base& a_rhs) const
    {
        bool res = shape() == a_rhs.shape();
        if (res)
            res = !memcmp(data(), a_rhs.data(), numel() * sizeof(T));
        return res;
    }

    template <typename TJS> operator TJS() const
    {
        return zax_to_json(0);
    }

    virtual std::string zax_to_json(int a_deep) const
    {
        unsigned int alloc_size = ZaxJsonParser::initial_alloc_size();
        char* json = new char[alloc_size];
        while (!zax_to_json(json, alloc_size - 1, a_deep))
            if (!ZaxJsonParser::reallocate_json(json, alloc_size))
                break;
        std::string a_json = json ? json : "";
        delete[] json;
        return a_json;
    }

    friend std::ostream& operator<<(std::ostream& os, const tensor_base& a_obj)
    {
        std::string s = a_obj;
        return os << s;
    }

    OPERATOR_EQ()
};

std::vector<int> get_dimensions(const char* a_json);

template<typename T>
inline void get_type_c(char* a_type_c);

template<>
inline void get_type_c<int>(char* a_type_c)
{
    strcpy(a_type_c, "%d,");
}

template<>
inline void get_type_c<unsigned int>(char* a_type_c)
{
    strcpy(a_type_c, "%u,");
}

template<>
inline void get_type_c<unsigned short int>(char* a_type_c)
{
    strcpy(a_type_c, "%hu,");
}

template<>
inline void get_type_c<short int>(char* a_type_c)
{
    strcpy(a_type_c, "%hd,");
}

template<>
inline void get_type_c<unsigned char>(char* a_type_c)
{
    strcpy(a_type_c, "%d,");
}

template<>
inline void get_type_c<int8_t>(char* a_type_c)
{
    strcpy(a_type_c, "%c,");
}

template<>
inline void get_type_c<float>(char* a_type_c)
{
    strcpy(a_type_c, "%f,");
}

template<>
inline void get_type_c<double>(char* a_type_c)
{
    strcpy(a_type_c, "%f,");
}

static inline void indent_new_line(char*& a_json, int& a_result, int a_deep)
{
    *(a_json + a_result++) = '\n';
    if (int nrindent = ZaxJsonParser::get_nr_indent())
        for (int i = 0; i < nrindent * a_deep; ++i)
            *(a_json + a_result++) = ' ';
    *(a_json + a_result++) = '[';
    *(a_json + a_result) = 0;
}

template<typename T>
struct tensor_1d: public tensor_base<T>
{
    using tensor_base<T>::d1;

public:
    T* m_data = 0;
    int m_wrap_around_bytes;

    tensor_1d(const tensor_1d& a_rhs)
    {
        *this = a_rhs;
    }

    tensor_1d(int a_size = 0, T* a_data = 0, int a_wrap_around_bytes = false)
        : m_wrap_around_bytes(a_wrap_around_bytes)
    {
        resize(a_size, a_data, a_wrap_around_bytes);
    }

    tensor_1d(const char* a_json)
    {
        zax_from_json(a_json);
    }

    tensor_1d(const std::string& a_json)
    {
        zax_from_json(a_json.c_str());
    }

    virtual ~tensor_1d()
    {
        destroy();
    }

    virtual tensor_1d& operator = (const tensor_1d& a_rhs)
    {
        resize(a_rhs.d1, a_rhs.data());
        return *this;
    }

    void destroy()
    {
        if (m_data && !m_wrap_around_bytes)
        {
            delete[] m_data;
            m_data = 0;
        }
    }

    inline const int size() const
    {
        return d1;
    }

    inline void clear()
    {
        resize(0);
    }

    virtual bool reshape(const std::vector<int>& a_sizes)
    {
        bool res = false;
        if (a_sizes.size() == 1 && d1 == a_sizes[0])
            res = true;
        return res;
    }

    virtual void resize(const std::vector<int>& a_sizes, T* a_data = 0, int a_wrap_around_bytes = false, bool a_reshaping = false)
    {
        if (a_sizes.size() > 0)
            resize(a_sizes[0], a_data, a_wrap_around_bytes);
    }

    void resize(int a_size, T* a_data = 0, int a_wrap_around_bytes = false)
    {
        if(a_size == d1 && !a_wrap_around_bytes)
        {
            if (a_data)
                memcpy(m_data, a_data, d1 * sizeof(T));
        }
        else
        {
            destroy();
            if (!a_size)
            {
                d1 = 0;
                m_data = 0;
            }
            else
            {
                d1 = a_size;
                if (a_data && a_wrap_around_bytes)
                    m_data = a_data;
                else
                {
                    m_data = new T[d1];
                    memset(m_data, 0, d1 * sizeof(T));
                }
                if (a_data && !a_wrap_around_bytes)
                    memcpy(m_data, a_data, d1 * sizeof(T));
            }
        }
        m_wrap_around_bytes = a_wrap_around_bytes;
    }

    virtual T* data() const
    {
        return m_data;
    }

    virtual void detach_data()
    {
        m_wrap_around_bytes = true;
    }

    virtual std::vector<int> shape() const
    {
        return {d1};
    }

    virtual void fill(const T& a_val)
    {
        for (int j = 0; j < d1; j++)
            m_data[j] = a_val;
    }

    T& operator [](int a_idx)
    {
        return m_data[a_idx];
    }

    T operator [](int a_idx) const
    {
        return m_data[a_idx];
    }

    virtual void zax_from_json(const char* a_json, std::string* a_err_stream = 0)
    {
        resize(get_dimensions(a_json));
        a_json = strchr(a_json, '[');
        if (a_json)
            for (int i = 0; i < d1; ++i)
            {
                char* endptr;
                m_data[i] = strtof(++a_json, &endptr);
                a_json = strchr(a_json, ',');
                if (!a_json)
                    break;
            }
    }

    virtual int zax_to_json(char* a_json, int a_alloc_size, int a_deep) const
    {
        char typed[8];
        get_type_c<T>(typed);

        a_json[0] = 0;
        strcat(a_json, "[");
        int offset = 1;
        for (int i = 0; i < d1; ++i)
            offset += snprintf(a_json + offset, a_alloc_size, typed, m_data[i]);
        sprintf(a_json + --offset, "%s", "]");
        if (!(offset))
        {
            sprintf(a_json, "[]");
            offset = 2;
        }
        return offset;
    }
};

template<typename T>
struct tensor_2d: public tensor_base<T>
{
    using tensor_base<T>::d1;
    using tensor_base<T>::d2;

public:
    T** m_data = 0;
    int m_wrap_around_bytes;

    tensor_2d(int a_d1 = 0, int a_d2 = 0, T* a_data = 0, int a_wrap_around_bytes = false)
        : m_wrap_around_bytes(a_wrap_around_bytes)
    {
        resize(a_d1, a_d2, a_data, a_wrap_around_bytes);
    }

    tensor_2d(const tensor_2d& a_rhs)
    {
        *this = a_rhs;
    }

    tensor_2d(const char* a_json)
    {
        zax_from_json(a_json);
    }

    tensor_2d(const std::string& a_json)
    {
        zax_from_json(a_json.c_str());
    }

    virtual ~tensor_2d()
    {
        destroy();
    }

    virtual tensor_2d& operator = (const tensor_2d& a_rhs)
    {
        resize(a_rhs.d1, a_rhs.d2, a_rhs.numel() ? a_rhs.data() : 0);
        return *this;
    }

    void destroy()
    {
        if (m_data)
        {
            if (!m_wrap_around_bytes)
                delete[] m_data[0];
            delete[] m_data;
            m_data = 0;
        }
    }

    bool reshape(int a_d1, int a_d2)
    {
        bool res = a_d1 * a_d2 == d1 * d2;
        if (res && (a_d1 != d1))
            resize(a_d1, a_d2, data(), false, true);
        return res;
    }

    virtual bool reshape(const std::vector<int>& a_sizes)
    {
        bool res = false;
        if (a_sizes.size() == 2)
            res = reshape(a_sizes[0], a_sizes[1]);
        return res;
    }

    virtual void resize(const std::vector<int>& a_sizes, T* a_data = 0, int a_wrap_around_bytes = false, bool a_reshaping = false)
    {
        if (a_sizes.size() > 1)
            resize(a_sizes[0], a_sizes[1], a_data, a_wrap_around_bytes, a_reshaping);
    }

    void resize(int a_d1, int a_d2, T* a_data = 0, int a_wrap_around_bytes = false, bool a_reshaping = false)
    {
        if(a_d1 == d1 && a_d2 == d2 && !a_wrap_around_bytes)
        {
            if (a_data)
                memcpy(&m_data[0][0], a_data, d1 * d2 * sizeof(T));
        }
        else
        {
            if (!a_reshaping)
                destroy();
            if ((!a_d1) || (!a_d2))
            {
                d2 = 0;
                d1 = 0;
                m_data = 0;
            }
            else
            {
                d2 = a_d2;
                d1 = a_d1;
                if (a_reshaping)
                    delete[] m_data;
                m_data = new T*[d1];
                if (a_data && (a_wrap_around_bytes || a_reshaping))
                    m_data[0] = a_data;
                else
                {
                    m_data[0] = new T[d1 * d2];
                    memset(m_data[0], 0, d1 * d2 * sizeof(T));
                }
                for (int i = 0; i < d1; ++i)
                    m_data[i] = &m_data[0][i * d2];
                if (a_data && !a_wrap_around_bytes && !a_reshaping)
                    memcpy(&m_data[0][0], a_data, d1 * d2 * sizeof(T));
            }
        }
        m_wrap_around_bytes = a_wrap_around_bytes;
    }

    virtual T* data() const
    {
        return m_data ? m_data[0] : 0;
    }

    virtual void detach_data()
    {
        m_wrap_around_bytes = true;
    }

    virtual std::vector<int> shape() const
    {
        return {d1, d2};
    }

    virtual void fill(const T& a_val)
    {
        for (int i = 0; i < d1; ++i)
            for (int j = 0; j < d2; ++j)
                m_data[i][j] = a_val;
    }

    T*& operator [](int a_idx)
    {
        return m_data[a_idx];
    }

    T* operator [](int a_idx) const
    {
        return m_data[a_idx];
    }

    virtual void zax_from_json(const char* a_json, std::string* a_err_stream = 0)
    {
        resize(get_dimensions(a_json));
        a_json = strchr(a_json, '[');
        if (a_json)
            for (int i = 0; i < d1; ++i)
                if ((a_json = strchr(++a_json, '[')))
                    for (int j = 0; j < d2; ++j)
                    {
                        char* endptr;
                        m_data[i][j] = strtof(++a_json, &endptr);
                        a_json = strchr(a_json, ',');
                        if (!a_json)
                            break;
                    }
    }

    virtual int zax_to_json(char* a_json, int a_alloc_size, int a_deep) const
    {
        char typed[8];
        get_type_c<T>(typed);

        a_json[0] = 0;
        int offset = 0;
        strcat(a_json, "[");
        offset++;
        for (int i = 0; i < d1; ++i)
        {
            indent_new_line(a_json, offset, a_deep + 1);
            for (int j = 0; j < d2; ++j)
                offset += snprintf(a_json + offset, a_alloc_size, typed, m_data[i][j]);
            offset += sprintf(a_json + --offset, "%s", "],");
        }
        sprintf(a_json + --offset, "%s", "]");
        if (!(--offset))
        {
            sprintf(a_json, "[[]]");
            offset = 4;
        }
        return offset;
    }
};

template<typename T>
struct tensor_3d: public tensor_base<T>
{
    using tensor_base<T>::d1;
    using tensor_base<T>::d2;
    using tensor_base<T>::d3;

public:
    T*** m_data = 0;
    int m_wrap_around_bytes;

    tensor_3d(const tensor_3d& a_rhs)
    {
        *this = a_rhs;
    }

    tensor_3d(int a_d1 = 0, int a_d2 = 0, int a_d3 = 0, T* a_data = 0, int a_wrap_around_bytes = false)
        : m_wrap_around_bytes(a_wrap_around_bytes)
    {
        resize(a_d1, a_d2, a_d3, a_data, a_wrap_around_bytes);
    }

    tensor_3d(const char* a_json)
    {
        zax_from_json(a_json);
    }

    tensor_3d(const std::string& a_json)
    {
        zax_from_json(a_json.c_str());
    }

    virtual ~tensor_3d()
    {
        destroy();
    }

    virtual tensor_3d& operator = (const tensor_3d& a_rhs)
    {
        resize(a_rhs.d1, a_rhs.d2, a_rhs.d3, ((tensor_3d*)&a_rhs)->numel() ? ((tensor_3d*)&a_rhs)->data() : 0);
        return *this;
    }

    void destroy()
    {
        if (m_data)
        {
            if (!m_wrap_around_bytes)
                delete[] m_data[0][0];
            delete[] m_data;
            m_data = 0;
        }
    }

    bool reshape(int a_d1, int a_d2, int a_d3)
    {
        bool res = a_d1 * a_d2 * a_d3 == d1 * d2 * d3;
        if (res && (a_d1 != d1 || a_d2 != d2))
            resize(a_d1, a_d2, a_d3, data(), false, true);
        return res;
    }

    virtual bool reshape(const std::vector<int>& a_sizes)
    {
        bool res = false;
        if (a_sizes.size() == 3)
            res = reshape(a_sizes[0], a_sizes[1], a_sizes[2]);
        return res;
    }

    virtual void resize(const std::vector<int>& a_sizes, T* a_data = 0, int a_wrap_around_bytes = false, bool a_reshaping = false)
    {
        if (a_sizes.size() > 2)
            resize(a_sizes[0], a_sizes[1], a_sizes[2], a_data, a_wrap_around_bytes, a_reshaping);
    }

    void resize(int a_d1, int a_d2, int a_d3, T* a_data = 0, int a_wrap_around_bytes = false, bool a_reshaping = false)
    {
        if(a_d1 == d1 && a_d2 == d2 && a_d3 == d3 && !a_wrap_around_bytes)
        {
            if (a_data)
                memcpy(m_data[0][0], a_data, d1 * d2 * d3 * sizeof(T));
        }
        else
        {
            if (!a_reshaping)
                destroy();
            if ((!a_d1) || (!a_d2) || (!a_d3))
            {
                d2 = 0;
                d1 = 0;
                d3 = 0;
                m_data = 0;
            }
            else
            {
                d2 = a_d2;
                d1 = a_d1;
                d3 = a_d3;

                if (a_reshaping)
                    delete[] m_data;
                m_data = new T**[d1 + d1 * d2];
                T** all_y = (T**)&m_data[d1];

                for (int x = 0; x < d1; ++x, all_y += d2)
                    m_data[x] = all_y;

                if (a_data && (a_wrap_around_bytes || a_reshaping))
                    m_data[0][0] = a_data;
                else
                {
                    m_data[0][0] = new T[d1 * d2 * d3];
                    memset(m_data[0][0], 0, d1 * d2 * d3 * sizeof(T));
                }

                if (a_data && !a_wrap_around_bytes && !a_reshaping)
                    memcpy(m_data[0][0], a_data, d1 * d2 * d3 * sizeof(T));

                T* all_r = m_data[0][0];
                for (int x = 0; x < d1; ++x, all_y += d2)
                    for (int y = 0; y < d2; ++y, all_r += d3)
                        m_data[x][y] = all_r;
            }
        }
        m_wrap_around_bytes = a_wrap_around_bytes;
    }

    virtual T* data() const
    {
        return m_data ? m_data[0][0] : 0;
    }

    virtual void detach_data()
    {
        m_wrap_around_bytes = true;
    }

    virtual std::vector<int> shape() const
    {
        return {d1, d2, d3};
    }

    virtual void fill(const T& a_val)
    {
        for (int i = 0; i < d1; i++)
            for (int j = 0; j < d2; j++)
                for (int k = 0; k < d3; k++)
                    m_data[i][j][k] = a_val;
    }

    T**& operator [](int a_idx)
    {
        return m_data[a_idx];
    }

    T** operator [](int a_idx) const
    {
        return m_data[a_idx];
    }

    virtual void zax_from_json(const char* a_json, std::string* a_err_stream = 0)
    {
        resize(get_dimensions(a_json));
        a_json = strchr(a_json, '[');
        if (a_json)
            for (int i = 0; i < d1; ++i)
                if ((a_json = strchr(++a_json, '[')))
                    for (int j = 0; j < d2; ++j)
                        if ((a_json = strchr(++a_json, '[')))
                            for (int k = 0; k < d3; ++k)
                            {
                                char* endptr;
                                m_data[i][j][k] = strtof(++a_json, &endptr);
                                a_json = strchr(a_json, ',');
                                if (!a_json)
                                    break;
                            }
    }

    virtual int zax_to_json(char* a_json, int a_alloc_size, int a_deep) const
    {
        char typed[8];
        get_type_c<T>(typed);

        a_json[0] = 0;
        int offset = 0;
        strcat(a_json, "[");
        offset++;
        for (int i = 0; i < d1; ++i)
        {
            indent_new_line(a_json, offset, a_deep + 1);
            for (int j = 0; j < d2; ++j)
            {
                indent_new_line(a_json, offset, a_deep + 2);
                for (int k = 0; k < d3; ++k)
                    offset += snprintf(a_json + offset, a_alloc_size, typed, m_data[i][j][k]);
                offset += sprintf(a_json + --offset, "%s", "],");
            }
            offset += sprintf(a_json + --offset, "%s", "],");
        }
        sprintf(a_json + --offset, "%s", "]");
        if (!(--offset))
        {
            sprintf(a_json, "[[[]]]");
            offset = 6;
        }
        return offset;
    }
};

template<typename T>
struct tensor_4d: public tensor_base<T>
{
    using tensor_base<T>::d1;
    using tensor_base<T>::d2;
    using tensor_base<T>::d3;
    using tensor_base<T>::d4;

public:
    T**** m_data = 0;
    int m_wrap_around_bytes;

    tensor_4d(const tensor_4d& a_rhs)
    {
        *this = a_rhs;
    }

    tensor_4d(int a_d1 = 0, int a_d2 = 0, int a_d3 = 0, int a_d4 = 0, T* a_data = 0, int a_wrap_around_bytes = false)
        : m_wrap_around_bytes(a_wrap_around_bytes)
    {
        resize(a_d1, a_d2, a_d3, a_d4, a_data, a_wrap_around_bytes);
    }

    tensor_4d(const char* a_json)
    {
        zax_from_json(a_json);
    }

    tensor_4d(const std::string& a_json)
    {
        zax_from_json(a_json.c_str());
    }

    virtual ~tensor_4d()
    {
        destroy();
    }

    virtual tensor_4d& operator = (const tensor_4d& a_rhs)
    {
        resize(a_rhs.d1, a_rhs.d2, a_rhs.d3, a_rhs.d4, ((tensor_4d*)&a_rhs)->numel() ? ((tensor_4d*)&a_rhs)->data() : 0);
        return *this;
    }

    void destroy()
    {
        if (m_data)
        {
            if (!m_wrap_around_bytes)
                delete[] m_data[0][0][0];
            delete[] m_data;
            m_data = 0;
        }
    }

    bool reshape(int a_d1, int a_d2, int a_d3, int a_d4)
    {
        bool res = a_d1 * a_d2 * a_d3 * a_d4 == d1 * d2 * d3 * d4;
        if (res && (a_d1 != d1 || a_d2 != d2 || a_d3 != d3))
            resize(a_d1, a_d2, a_d3, a_d4, data(), false, true);
        return res;
    }

    virtual bool reshape(const std::vector<int>& a_sizes)
    {
        bool res = false;
        if (a_sizes.size() == 4)
            res = reshape(a_sizes[0], a_sizes[1], a_sizes[2], a_sizes[3]);
        return res;
    }

    virtual void resize(const std::vector<int>& a_sizes, T* a_data = 0, int a_wrap_around_bytes = false, bool a_reshaping = false)
    {
        if (a_sizes.size() > 3)
            resize(a_sizes[0], a_sizes[1], a_sizes[2], a_sizes[3], a_data, a_wrap_around_bytes, a_reshaping);
    }

    void resize(int a_d1, int a_d2, int a_d3, int a_d4, T* a_data = 0, int a_wrap_around_bytes = false, bool a_reshaping = false)
    {
        if(a_d1 == d1 && a_d2 == d2 && a_d3 == d3 && a_d4 == d4 && !a_wrap_around_bytes)
        {
            if (a_data)
                memcpy(m_data[0][0][0], a_data, d1 * d2 * d3 * d4 * sizeof(T));
        }
        else
        {
            if (!a_reshaping)
                destroy();
            if ((!a_d1) || (!a_d2) || (!a_d3) || (!a_d4))
            {
                d2 = 0;
                d1 = 0;
                d4 = 0;
                d3 = 0;
                m_data = 0;
            }
            else
            {
                d2 = a_d2;
                d1 = a_d1;
                d4 = a_d4;
                d3 = a_d3;

                if (a_reshaping)
                    delete[] m_data;
                m_data = new T***[d1 + d1 * d2 + d1 * d2 * d3];
                T*** all_y = (T***)&m_data[d1];
                T** all_r = (T**)&m_data[d1 + d1 * d2];

                for (int x = 0; x < d1; ++x, all_y += d2)
                {
                    m_data[x] = all_y;
                    for (int y = 0; y < d2; ++y, all_r += d3)
                        m_data[x][y] = all_r;
                }

                if (a_data && (a_wrap_around_bytes || a_reshaping))
                    m_data[0][0][0] = a_data;
                else
                {
                    m_data[0][0][0] = new T[d1 * d2 * d3 * d4];
                    memset(m_data[0][0][0], 0, d1 * d2 * d3 * d4 * sizeof(T));
                }

                if (a_data && !a_wrap_around_bytes && !a_reshaping)
                    memcpy(m_data[0][0][0], a_data, d1 * d2 * d3 * d4 * sizeof(T));

                T* all_c = m_data[0][0][0];
                for (int x = 0; x < d1; ++x, all_y += d2)
                    for (int y = 0; y < d2; ++y, all_r += d3)
                        for (int r = 0; r < d3; ++r, all_c += d4)
                            m_data[x][y][r] = all_c;
            }
        }
        m_wrap_around_bytes = a_wrap_around_bytes;
    }

    virtual T* data() const
    {
        return m_data ? m_data[0][0][0] : 0;
    }

    virtual void detach_data()
    {
        m_wrap_around_bytes = true;
    }

    virtual std::vector<int> shape() const
    {
        return {d1, d2, d3, d4};
    }

    virtual void fill(const T& a_val)
    {
        for (int i = 0; i < d1; i++)
            for (int j = 0; j < d2; j++)
                for (int k = 0; k < d3; k++)
                    for (int l = 0; l < d4; l++)
                        m_data[i][j][k][l] = a_val;
    }

    T***& operator [](int a_idx)
    {
        return m_data[a_idx];
    }

    T*** operator [](int a_idx) const
    {
        return m_data[a_idx];
    }

    virtual void zax_from_json(const char* a_json, std::string* a_err_stream = 0)
    {
        resize(get_dimensions(a_json));
        a_json = strchr(a_json, '[');
        if (a_json)
            for (int i = 0; i < d1; ++i)
                if ((a_json = strchr(++a_json, '[')))
                    for (int j = 0; j < d2; ++j)
                        if ((a_json = strchr(++a_json, '[')))
                            for (int k = 0; k < d3; ++k)
                                if ((a_json = strchr(++a_json, '[')))
                                    for (int l = 0; l < d4; ++l)
                                    {
                                        char* endptr;
                                        m_data[i][j][k][l] = strtof(++a_json, &endptr);
                                        a_json = strchr(a_json, ',');
                                        if (!a_json)
                                            break;
                                    }
    }

    virtual int zax_to_json(char* a_json, int a_alloc_size, int a_deep) const
    {
        char typed[8];
        get_type_c<T>(typed);

        a_json[0] = 0;
        int offset = 0;
        strcat(a_json, "[");
        offset++;
        for (int i = 0; i < d1; ++i)
        {
            indent_new_line(a_json, offset, a_deep + 1);
            for (int j = 0; j < d2; ++j)
            {
                indent_new_line(a_json, offset, a_deep + 2);
                for (int k = 0; k < d3; ++k)
                {
                    indent_new_line(a_json, offset, a_deep + 3);
                    for (int l = 0; l < d4; ++l)
                        offset += snprintf(a_json + offset, a_alloc_size, typed, m_data[i][j][k][l]);
                    offset += sprintf(a_json + --offset, "%s", "],");
                }
                offset += sprintf(a_json + --offset, "%s", "],");
            }
            offset += sprintf(a_json + --offset, "%s", "],");
        }
        sprintf(a_json + --offset, "%s", "]");
        if (!(--offset))
        {
            sprintf(a_json, "[[[[]]]]");
            offset = 8;
        }
        return offset;
    }
};

class float_1d: public tensor_1d<float>
{
    using tensor_1d<float>::tensor_1d;
    OPERATOR_EQ()
};

class float_2d: public tensor_2d<float>
{
    using tensor_2d<float>::tensor_2d;
    OPERATOR_EQ()
};

class float_3d: public tensor_3d<float>
{
    using tensor_3d<float>::tensor_3d;
    OPERATOR_EQ()
};

class float_4d: public tensor_4d<float>
{
    using tensor_4d<float>::tensor_4d;
    OPERATOR_EQ()
};

class double_1d: public tensor_1d<double>
{
    using tensor_1d<double>::tensor_1d;
    OPERATOR_EQ()
};

class double_2d: public tensor_2d<double>
{
    using tensor_2d<double>::tensor_2d;
    OPERATOR_EQ()
};

class double_3d: public tensor_3d<double>
{
    using tensor_3d<double>::tensor_3d;
    OPERATOR_EQ()
};

class double_4d: public tensor_4d<double>
{
    using tensor_4d<double>::tensor_4d;
    OPERATOR_EQ()
};

class int_1d: public tensor_1d<int>
{
    using tensor_1d<int>::tensor_1d;
    OPERATOR_EQ()
};

class int_2d: public tensor_2d<int>
{
    using tensor_2d<int>::tensor_2d;
    OPERATOR_EQ()
};

class int_3d: public tensor_3d<int>
{
    using tensor_3d<int>::tensor_3d;
    OPERATOR_EQ()
};

class int_4d: public tensor_4d<int>
{
    using tensor_4d<int>::tensor_4d;
    OPERATOR_EQ()
};

class uint_1d: public tensor_1d<unsigned int>
{
    using tensor_1d<unsigned int>::tensor_1d;
    OPERATOR_EQ()
};

class uint_2d: public tensor_2d<unsigned int>
{
    using tensor_2d<unsigned int>::tensor_2d;
    OPERATOR_EQ()
};

class uint_3d: public tensor_3d<unsigned int>
{
    using tensor_3d<unsigned int>::tensor_3d;
    OPERATOR_EQ()
};

class uint_4d: public tensor_4d<unsigned int>
{
    using tensor_4d<unsigned int>::tensor_4d;
    OPERATOR_EQ()
};

class ui8_1d: public tensor_1d<unsigned char>
{
    using tensor_1d<unsigned char>::tensor_1d;
    OPERATOR_EQ()
};

class ui8_2d: public tensor_2d<unsigned char>
{
    using tensor_2d<unsigned char>::tensor_2d;
    OPERATOR_EQ()
};

class ui8_3d: public tensor_3d<unsigned char>
{
    using tensor_3d<unsigned char>::tensor_3d;
    OPERATOR_EQ()
};

class ui8_4d: public tensor_4d<unsigned char>
{
    using tensor_4d<unsigned char>::tensor_4d;
    OPERATOR_EQ()
};

class i8_1d: public tensor_1d<int8_t>
{
    using tensor_1d<int8_t>::tensor_1d;
    OPERATOR_EQ()
};

class i8_2d: public tensor_2d<int8_t>
{
    using tensor_2d<int8_t>::tensor_2d;
    OPERATOR_EQ()
};

class i8_3d: public tensor_3d<int8_t>
{
    using tensor_3d<int8_t>::tensor_3d;
    OPERATOR_EQ()
};

class i8_4d: public tensor_4d<int8_t>
{
    using tensor_4d<int8_t>::tensor_4d;
    OPERATOR_EQ()
};

class ui16_1d: public tensor_1d<unsigned short>
{
    using tensor_1d<unsigned short>::tensor_1d;
    OPERATOR_EQ()
};

class i16_1d: public tensor_1d<short>
{
    using tensor_1d<short>::tensor_1d;
    OPERATOR_EQ()
};

class ui16_2d: public tensor_2d<unsigned short>
{
    using tensor_2d<unsigned short>::tensor_2d;
    OPERATOR_EQ()
};

class i16_2d: public tensor_2d<short>
{
    using tensor_2d<short>::tensor_2d;
    OPERATOR_EQ()
};

class ui16_3d: public tensor_3d<unsigned short>
{
    using tensor_3d<unsigned short>::tensor_3d;
    OPERATOR_EQ()
};

class i16_3d: public tensor_3d<short>
{
    using tensor_3d<short>::tensor_3d;
    OPERATOR_EQ()
};

class ui16_4d: public tensor_4d<unsigned short>
{
    using tensor_4d<unsigned short>::tensor_4d;
    OPERATOR_EQ()
};

class i16_4d: public tensor_4d<short>
{
    using tensor_4d<short>::tensor_4d;
    OPERATOR_EQ()
};

template<typename T, typename T_1d, typename T_2d, typename T_3d, typename T_4d>
class tensor_t: public tensor_base<T>
{
    tensor_base<T>* current = 0;
public:
    using tensor_base<T>::numel;
    T_1d* m_1d = 0;
    T_2d* m_2d = 0;
    T_3d* m_3d = 0;
    T_4d* m_4d = 0;

    tensor_t()
    {}

    tensor_t(const char* a_json)
    {
        zax_from_json(a_json);
    }

    tensor_t(std::string& a_json)
    {
        zax_from_json(a_json.c_str());
    }

    tensor_t(const tensor_t& a_rhs)
    {
        resize(a_rhs.shape(), a_rhs.data());
    }

    tensor_t(const std::vector<int>& a_sizes, T* a_data = 0, int a_wrap_around_bytes = false)
    {
        resize(a_sizes, a_data, a_wrap_around_bytes);
    }

    virtual ~tensor_t()
    {
        delete current;
    }

    void create(const std::vector<int>& a_sizes, T* a_data = 0, int a_wrap_around_bytes = false)
    {
        if (current)
            delete current;
        current = 0;
        m_1d = 0;
        m_2d = 0;
        m_3d = 0;
        m_4d = 0;
        if (a_sizes.size() == 1)
            current = m_1d = new T_1d;
        else if (a_sizes.size() == 2)
            current = m_2d = new T_2d;
        else if (a_sizes.size() == 3)
            current = m_3d = new T_3d;
        else if (a_sizes.size() == 4)
            current = m_4d = new T_4d;
        if (current)
            current->resize(a_sizes, a_data, a_wrap_around_bytes);
        switch (a_sizes.size())
        {
        case 4:
            tensor_base<T>::d4 = a_sizes[3];
        case 3:
            tensor_base<T>::d3 = a_sizes[2];
        case 2:
            tensor_base<T>::d2 = a_sizes[1];
        case 1:
            tensor_base<T>::d1 = a_sizes[0];
        default:
            break;
        }
    }

    virtual T* data() const
    {
        if (!current)
            return 0;
        else
            return current->data();
    }

    virtual void fill(const T& a_val)
    {
        current->fill(a_val);
    }

    virtual T** data_2d() const
    {
        return m_2d->m_data;
    }

    virtual T*** data_3d() const
    {
        return m_3d->m_data;
    }

    virtual T**** data_4d() const
    {
        return m_4d->m_data;
    }

    virtual void detach_data()
    {
        if (current)
            current->detach_data();
    }

    virtual std::vector<int> shape() const
    {
        if (current)
            return current->shape();
        else
            return {};
    }

    void reshape(int d1, int d2 = -1, int d3 = -1, int d4 = -1)
    {
        std::vector<int> sizes = {d1};
        if (d2 != -1)
            sizes.push_back(d2);
        if (d3 != -1)
            sizes.push_back(d3);
        if (d4 != -1)
            sizes.push_back(d4);
        reshape(sizes);
    }

    virtual bool reshape(const std::vector<int>& a_sizes)
    {
        bool res = numel() == std::accumulate(a_sizes.begin(), a_sizes.end(), 1, std::multiplies<int>());
        if (res)
            resize(a_sizes, data(), false, true);
        return res;
    }

    void resize(int d1, int d2 = -1, int d3 = -1, int d4 = -1)
    {
        std::vector<int> sizes = {d1};
        if (d2 != -1)
            sizes.push_back(d2);
        if (d3 != -1)
            sizes.push_back(d3);
        if (d4 != -1)
            sizes.push_back(d4);
        resize(sizes);
    }

    virtual void resize(const std::vector<int>& a_sizes, T* a_data = 0, int a_wrap_around_bytes = false, bool a_reshaping = false)
    {
        switch (a_sizes.size())
        {
        case 4:
            tensor_base<T>::d4 = a_sizes[3];
        case 3:
            tensor_base<T>::d3 = a_sizes[2];
        case 2:
            tensor_base<T>::d2 = a_sizes[1];
        case 1:
            tensor_base<T>::d1 = a_sizes[0];
        default:
            break;
        }
        if (current)
        {
            if (numel() == std::accumulate(a_sizes.begin(), a_sizes.end(), 1, std::multiplies<int>()))
            {
                if (shape().size() == a_sizes.size())
                    current->resize(a_sizes, a_data ? a_data : data(), a_wrap_around_bytes, true);
                else
                {
                    T* dt = a_data ? a_data : data();
                    detach_data();
                    std::vector<int> ones(a_sizes.size());
                    std::fill(ones.begin(), ones.end(), 1);
                    create(ones);
                    current->resize(a_sizes, dt, true, true);
                }
            }
            else
                create(a_sizes, a_data, a_wrap_around_bytes);
        }
        else
            create(a_sizes, a_data, a_wrap_around_bytes);
        tensor_base<T>::d1d = data();
        if (m_2d)
            tensor_base<T>::d2d = data_2d();
        else
            tensor_base<T>::d2d = 0;
        if (m_3d)
            tensor_base<T>::d3d = data_3d();
        else
            tensor_base<T>::d3d = 0;
        if (m_4d)
            tensor_base<T>::d4d = data_4d();
        else
            tensor_base<T>::d4d = 0;
    }

    bool view(const std::vector<int>& a_sizes, tensor_t& a_orig)
    {
        bool res = a_orig.numel() == std::accumulate(a_sizes.begin(), a_sizes.end(), 1, std::multiplies<int>());
        if (res)
            resize(a_sizes, a_orig.data(), true);
        return res;
    }

    tensor_t& squeeze(unsigned int a_index_of_dimension, bool* success = 0)
    {
        if (success)
            *success = false;
        std::vector<int> sh = shape();
        if (sh.size() > a_index_of_dimension)
            if (sh[a_index_of_dimension] == 1)
            {
                sh.erase(sh.begin() + a_index_of_dimension);
                resize(sh, data(), false, true);
                if (success)
                    *success = true;
            }
        return *this;
    }

    tensor_t& unsqueeze(unsigned int a_index_of_dimension, bool* success = 0)
    {
        if (success)
            *success = false;
        std::vector<int> sh = shape();
        if (sh.size() >= a_index_of_dimension)
        {
            sh.insert(sh.begin() + a_index_of_dimension, 1);
            resize(sh, data(), false, true);
            if (success)
                *success = true;
        }
        return *this;
    }

    T& operator()(int d1)
    {
        return (*m_1d)[d1];
    }

    T& operator()(int d1, int d2)
    {
        return (*m_2d)[d1][d2];
    }

    T& operator()(int d1, int d2, int d3)
    {
        return (*m_3d)[d1][d2][d3];
    }

    T& operator()(int d1, int d2, int d3, int d4)
    {
        return (*m_4d)[d1][d2][d3][d4];
    }

    virtual tensor_t& operator = (const T_1d& a_rhs)
    {
        create(a_rhs.shape(), a_rhs.data());
        return *this;
    }
    virtual tensor_t& operator = (const T_2d& a_rhs)
    {
        create(a_rhs.shape(), a_rhs.data());
        return *this;
    }
    virtual tensor_t& operator = (const T_3d& a_rhs)
    {
        create(a_rhs.shape(), a_rhs.data());
        return *this;
    }
    virtual tensor_t& operator = (const T_4d& a_rhs)
    {
        create(a_rhs.shape(), a_rhs.data());
        return *this;
    }

    virtual void zax_from_json(const char* a_json, std::string* a_err_stream = 0)
    {
        std::vector<int> sizes = get_dimensions(a_json);
        if (shape().size() != sizes.size())
            create(sizes);
        current->zax_from_json(a_json, a_err_stream);
    }

    virtual int zax_to_json(char* a_json, int a_alloc_size, int a_deep) const
    {
        if (current)
            return current->zax_to_json(a_json, a_alloc_size, a_deep);
        else
        {
            strcpy(a_json, "[]");
            return 2;
        }
    }

    OPERATOR_EQ()
};

#define tensor_f32 tensor_t<float, float_1d, float_2d, float_3d, float_4d>
#define tensor_f64 tensor_t<double, double_1d, double_2d, double_3d, double_4d>
#define tensor_i32 tensor_t<int, int_1d, int_2d, int_3d, int_4d>
#define tensor_ui32 tensor_t<unsigned int, uint_1d, uint_2d, uint_3d, uint_4d>
#define tensor_ui8 tensor_t<unsigned char, ui8_1d, ui8_2d, ui8_3d, ui8_4d>
#define tensor_i8 tensor_t<int8_t, i8_1d, i8_2d, i8_3d, i8_4d>
#define tensor_ui16 tensor_t<unsigned short, ui16_1d, ui16_2d, ui16_3d, ui16_4d>
#define tensor_i16 tensor_t<short, i16_1d, i16_2d, i16_3d, i16_4d>

template<typename arr_t, typename arr_t1 = int, typename arr_t2 = int, typename arr_t3 = int>
class array_of_tensors
{
public:
    std::vector<arr_t> m_data;

    arr_t& operator [](int a_idx)
    {
        return m_data[a_idx];
    }

    arr_t operator [](int a_idx) const
    {
        return m_data[a_idx];
    }

    void resize(int a_size)
    {
        m_data.resize(a_size);
    }

    int size() const
    {
        return m_data.size();
    }

    bool operator == (const array_of_tensors& a_rhs) const
    {
        bool res = size() == a_rhs.size();
        if (res)
            for (unsigned int i = 0; i < m_data.size(); ++i)
                if (!(res = m_data[i] == a_rhs[i]))
                    break;
        return res;
    }

    ZAX_JSON_SERIALIZABLE(array_of_tensors, JSON_PROPERTY(m_data, "^"))
};

class array_of_tensor_f32: public array_of_tensors<tensor_f32>
{
    using array_of_tensors<tensor_f32>::array_of_tensors;
    OPERATOR_EQ()
};

class array_of_tensor_i32: public array_of_tensors<tensor_i32>
{
    using array_of_tensors<tensor_i32>::array_of_tensors;
    OPERATOR_EQ()
};

class array_of_tensor_ui32: public array_of_tensors<tensor_ui32>
{
    using array_of_tensors<tensor_ui32>::array_of_tensors;
    OPERATOR_EQ()
};

class array_of_tensor_ui8: public array_of_tensors<tensor_ui8>
{
    using array_of_tensors<tensor_ui8>::array_of_tensors;
    OPERATOR_EQ()
};

class array_of_tensor_i8: public array_of_tensors<tensor_i8>
{
    using array_of_tensors<tensor_i8>::array_of_tensors;
    OPERATOR_EQ()
};

class array_of_tensor_ui16: public array_of_tensors<tensor_ui16>
{
    using array_of_tensors<tensor_ui16>::array_of_tensors;
    OPERATOR_EQ()
};

class array_of_tensor_i16: public array_of_tensors<tensor_i16>
{
    using array_of_tensors<tensor_i16>::array_of_tensors;
    OPERATOR_EQ()
};

#endif /// _ZAX_TENSOR_H_
