#pragma once

#include "common.h"
#include "matrix.h"

using namespace ParallelEM;

class TProcessor {
public:
  TProcessor(std::shared_ptr<TFilenameQueue> filename_queue,
    std::shared_ptr<TPhiMatrix> phi,
    std::shared_ptr<TPhiMatrix> n_wt,
    std::shared_ptr<TPhiMatrix> n_t,
    float* peprplexity_value,
    int num_inner_iter);

  TProcessor(const TProcessor&) = delete;

  ~TProcessor();

  bool is_stopped() const { return is_stopped_; }

private:
  std::shared_ptr<TFilenameQueue> filename_queue_;
  std::shared_ptr<TPhiMatrix> phi_;
  std::shared_ptr<TPhiMatrix> n_wt_;
  std::shared_ptr<TPhiMatrix> n_t_;
  float* perplexity_value_;
  int num_inner_iter_;

  mutable std::atomic<bool> is_stopped_;
  boost::thread thread_;

  void ThreadFunction();
};
