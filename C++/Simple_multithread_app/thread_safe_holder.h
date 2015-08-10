
#ifndef WC_THREAD_SAFE_HOLDER_H_
#define WC_THREAD_SAFE_HOLDER_H_

#include <unordered_map>
#include <memory>
#include <vector>
#include <queue>

#include "boost/thread/locks.hpp"
#include "boost/thread/mutex.hpp"


namespace word_counter_ns {

template<typename T>
class ThreadSafeHolder : boost::noncopyable {
 public:
  explicit ThreadSafeHolder()
      : lock_(), object_(std::make_shared<T>()) {}

  explicit ThreadSafeHolder(const std::shared_ptr<T>& object)
      : lock_(), object_(object) {}

  ~ThreadSafeHolder() {}

  std::shared_ptr<T> get() const {
    boost::lock_guard<boost::mutex> guard(lock_);
    return object_;
  }

  std::shared_ptr<T> get_copy() const {
    boost::lock_guard<boost::mutex> guard(lock_);
    return std::make_shared<T>(*object_);
  }

  void set(const std::shared_ptr<T>& object) {
    boost::lock_guard<boost::mutex> guard(lock_);
    object_ = object;
  }

 private:
  mutable boost::mutex lock_;
  std::shared_ptr<T> object_;
};

template<typename K, typename T>
class ThreadSafeCollectionHolder : boost::noncopyable {
 public:
  std::shared_ptr<T> get(const K& key) const {
    boost::lock_guard<boost::mutex> guard(lock_);
    return get_locked(key);
  }

  bool has_key(const K& key) const {
    boost::lock_guard<boost::mutex> guard(lock_);
    return object_.find(key) != object_.end();
  }

  void erase(const K& key) {
    boost::lock_guard<boost::mutex> guard(lock_);
    auto iter = object_.find(key);
    if (iter != object_.end()) {
      object_.erase(iter);
    }
  }

  std::shared_ptr<T> get_copy(const K& key) const {
    boost::lock_guard<boost::mutex> guard(lock_);
    auto value = get_locked(key);
    return value != nullptr ? std::make_shared<T>(*value) : std::shared_ptr<T>();
  }

  void set(const K& key, const std::shared_ptr<T>& value) {
    boost::lock_guard<boost::mutex> guard(lock_);
    auto iter = object_.find(key);
    if (iter == object_.end()) {
      object_.insert(std::pair<K, std::shared_ptr<T> >(key, value));
    } else {
      iter->second = value;
    }
  }

  void clear() {
    boost::lock_guard<boost::mutex> guard(lock_);
    object_.clear();
  }

  std::vector<K> keys() const {
    boost::lock_guard<boost::mutex> guard(lock_);
    std::vector<K> keys;
    for (auto iter = object_.begin(); iter != object_.end(); ++iter) {
      keys.push_back(iter->first);
    }
    return keys;
  }

 private:
  std::unordered_map<K, std::shared_ptr<T>> object_;
  mutable boost::mutex lock_;

  std::shared_ptr<T> get_locked(const K& key) const {
    auto iter = object_.find(key);
    return (iter == object_.end()) ? std::shared_ptr<T>() : iter->second;
  }
};

template<typename T>
class ThreadSafeQueue : boost::noncopyable {
 public:
  bool try_pop(T* elem) {
    boost::lock_guard<boost::mutex> guard(lock_);
    if (queue_.empty())
      return false;

    T tmp_elem = queue_.front();
    queue_.pop();
    *elem = tmp_elem;
    return true;
  }

  void push(const T& elem) {
    boost::lock_guard<boost::mutex> guard(lock_);
    queue_.push(elem);
  }

  int size() const {
    boost::lock_guard<boost::mutex> guard(lock_);
    return queue_.size();
  }

  int empty() const {
    boost::lock_guard<boost::mutex> guard(lock_);
    return queue_.empty();
  }

 private:
  mutable boost::mutex lock_;
  std::queue<T> queue_;
};

} // namespace word_counter_ns

#endif // WC_THREAD_SAFE_HOLDER_H_