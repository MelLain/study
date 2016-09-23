#pragma once

#include "common.h"

namespace ParallelEM {

  class TPhiMatrix {
  public:
    TPhiMatrix(int num_topics, const std::vector<Token>& tokens, bool random = true);

    TPhiMatrix(const TPhiMatrix&) = delete;

    float get_value(const Token& token, int topic_id) const;
    void add_value(const Token& token, int topic_id, float value);

    float get_value_unsafe(const Token& token, int topic_id) const;
    void set_value_unsafe(const Token& token, int topic_id, float value);

    int num_tokens() const { return num_tokens_; }
    int num_topics() const { return num_topics_; }

  private:
    int num_topics_;
    int num_tokens_;
    std::unordered_map<Token, int> token2id_;
    std::vector<std::vector<float> > data_;
    std::vector<std::shared_ptr<boost::mutex> > lockers_;
  };

}  // namespace ParallelEM
