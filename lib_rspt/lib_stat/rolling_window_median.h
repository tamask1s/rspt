/***************************************************************************
* Copyright 2024 Tamas Levente Kis - tamkis@gmail.com                      *
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

/**
 * Rolling median calculator using std::multiset as a heap.
 * Maintains a fixed-size sliding window and calculates the median
 * efficiently by moving an iterator rather than recalculating it every time.
 * The rolling_window_median class computes the median of a rolling window
 * of fixed size over a stream of data. A "rolling window" maintains
 * the last 'size_' elements and discards the oldest one as new elements are added.
 * The median shall stay correct even for duplicate values.
 */

template <typename T>
class preallocated_allocator
{
    T* mem_pool_ = nullptr;
    size_t nr_used_ = 0;
    std::vector<T*> free_addresses_;

public:
    using value_type = T;
    size_t capacity = 0;

    preallocated_allocator() = default;

    explicit preallocated_allocator(size_t max_elements)
        :nr_used_(0),
         capacity(max_elements)
    {
        if (capacity > 0)
            mem_pool_ = static_cast<T*>(::operator new(capacity * sizeof(T)));
    }

    preallocated_allocator(const preallocated_allocator& other)
        :nr_used_(other.nr_used_),
         free_addresses_(other.free_addresses_),
         capacity(other.capacity)
    {
        if (capacity > 0)
        {
            mem_pool_ = static_cast<T*>(::operator new(capacity * sizeof(T)));
            std::copy(other.mem_pool_, other.mem_pool_ + nr_used_, mem_pool_);
        }
    }

    preallocated_allocator(preallocated_allocator&& other) noexcept
        :mem_pool_(other.mem_pool_),
         nr_used_(other.nr_used_),
         free_addresses_(std::move(other.free_addresses_)),
         capacity(other.capacity)
    {
        other.mem_pool_ = nullptr;
        other.capacity = 0;
        other.nr_used_ = 0;
    }

    template <typename U>
    preallocated_allocator(const preallocated_allocator<U>& other) noexcept
        :nr_used_(0),
         capacity(other.capacity)
    {
        if (capacity > 0)
            mem_pool_ = static_cast<T*>(::operator new(capacity * sizeof(T)));
    }

    template <typename U>
    preallocated_allocator(preallocated_allocator<U>&& other) noexcept
        :nr_used_(0),
         capacity(other.capacity)
    {
        if (capacity > 0)
            mem_pool_ = static_cast<T*>(::operator new(capacity * sizeof(T)));
    }

    ~preallocated_allocator()
    {
        if (mem_pool_)
            ::operator delete(mem_pool_);
    }

    T* allocate(size_t n)
    {
        if (capacity == 0 || n > capacity)
            throw std::bad_alloc();

        if (free_addresses_.size() >= n)
        {
            T* ptr = free_addresses_.back();
            free_addresses_.pop_back();
            return ptr;
        }

        if (nr_used_ + n > capacity)
            throw std::bad_alloc();

        T* ptr = &mem_pool_[nr_used_];
        nr_used_ += n;
        return ptr;
    }

    void deallocate(T* ptr, size_t n)
    {
        for (size_t i = 0; i < n; ++i)
            free_addresses_.push_back(ptr + i);
    }

    template <typename U, typename... Args>
    void construct(U* ptr, Args&&... args)
    {
        new (ptr) U(std::forward<Args>(args)...);
    }

    template <typename U>
    void destroy(U* ptr)
    {
        ptr->~U();
    }

    template <typename U>
    struct rebind
    {
        using other = preallocated_allocator<U>;
    };

    bool operator==(const preallocated_allocator& other) const noexcept
    {
        return mem_pool_ == other.mem_pool_ && capacity == other.capacity;
    }

    bool operator!=(const preallocated_allocator& other) const noexcept
    {
        return !(*this == other);
    }
};

#define ODD (heap_.size() % 2)
template<typename T>
class rolling_window_median
{
    using Alloc = preallocated_allocator<T>;
    size_t size_;
    std::deque<typename std::multiset<T>::iterator> iterator_ring_;
    std::multiset<T, std::less<T>, Alloc> heap_;
    typename std::multiset<T>::iterator median_it_;

public:
    rolling_window_median(size_t size)
        :size_(size),
         heap_(std::less<T>(), Alloc(size + 1))
    {}

    /**
     * Inserts a new value into the sliding window,
     * and removes the oldest value (if necessary).
     * Updates the median iterator.
     * @param value: The new value to be added to the sliding window.
     * @return: The current median after inserting the value.
     */
    T insert(T value)
    {
        iterator_ring_.push_back(heap_.insert(value));

        if (iterator_ring_.size() > size_)
        {
            const auto& oldest_it_ = iterator_ring_.front();
            const T median_value_ = *median_it_; // copy is more optimal for bigger window sizes
            const T oldest_value_ = *oldest_it_;

            if (value < median_value_)
            {
                if (oldest_value_ < median_value_)
                {
                    heap_.erase(oldest_it_);
                }
                else if (ODD)
                {
                    if (oldest_value_ == median_value_)
                        ++median_it_;
                    heap_.erase(oldest_it_);
                    if (!ODD)
                        --median_it_;
                }
                else
                {
                    --median_it_;
                    if (oldest_value_ == *median_it_)
                        ++median_it_;
                    heap_.erase(oldest_it_);
                    if (!ODD)
                        --median_it_;
                }
            }
            else
            {
                if ((value > median_value_) && (oldest_value_ > median_value_))
                {
                    heap_.erase(oldest_it_);
                }
                else if (ODD)
                {
                    ++median_it_;
                    if (oldest_value_ <= *median_it_)
                        ++median_it_;
                    heap_.erase(oldest_it_);
                    if (!ODD)
                        --median_it_;
                }
                else
                {
                    if (oldest_value_ <= median_value_)
                        ++median_it_;
                    heap_.erase(oldest_it_);
                    if (!ODD)
                        --median_it_;
                }
            }
            iterator_ring_.pop_front();
        }
        else if (heap_.size() > 1)
        {
            if (value < *median_it_)
            {
                if (!ODD)
                    --median_it_;
            }
            else
            {
                if (ODD)
                    ++median_it_;
            }
        }
        else
            median_it_ = heap_.begin();

        if (!ODD)
            return (*median_it_ + *next(median_it_)) / 2.0;
        else
            return *median_it_;
    }
};
