#include "matrix.h"

namespace ParallelEM {

  TPhiMatrix::TPhiMatrix(int num_topics, const std::vector<std::string>& tokens, bool random)
      : num_topics_(num_topics)
      , num_tokens_(tokens.size()) {
    std::vector<float> normalizer(num_topics, 0.0f);

    // initialization
    for (int token_id = 0; token_id < tokens.size(); ++token_id) {
      token2id_.insert(std::make_pair(tokens[token_id], token_id));
      std::vector<float> token_row;
      for (int topic_id = 0; topic_id < num_topics_; ++topic_id) {
        float value = random ? rand() : 0.0f;
        token_row.push_back(value);
        normalizer[topic_id] += value;
      }
      lockers_.push_back(std::shared_ptr<boost::mutex>(new boost::mutex()));
      data_.push_back(std::move(token_row));
    }

    // normalization
    for (int token_id = 0; token_id < token2id_.size(); ++token_id) {
      for (int topic_id = 0; topic_id < num_topics; ++topic_id) {
        data_[token_id][topic_id] = normalizer[topic_id] > 0.0f ? data_[token_id][topic_id] / normalizer[topic_id] : 0.0f;
      }
    }
  }

  float TPhiMatrix::get_value(const Token& token, int topic_id) const {
    auto iter = token2id_.find(token);

    if (iter == token2id_.end() || topic_id >= num_topics_) {
      throw std::runtime_error(std::string("TPhiMatrix::get_value: index outside") +
        " bounds or no such token in matrix. Token: " + token + ", topic_id: " + std::to_string(topic_id));
    }

    int token_id = iter->second;
    boost::lock_guard<boost::mutex> guard(*lockers_[token_id]);
    return data_[token_id][topic_id];
  }

  void TPhiMatrix::add_value(const Token& token, int topic_id, float value) {
    auto iter = token2id_.find(token);

    if (iter == token2id_.end() || topic_id >= num_topics_) {
      throw std::runtime_error(std::string("TPhiMatrix::set_value: index outside") +
        " bounds or no such token in matrix. Token: " + token + ", topic_id: " + std::to_string(topic_id));
    }

    int token_id = iter->second;
    boost::lock_guard<boost::mutex> guard(*lockers_[token_id]);
    data_[token_id][topic_id] += value;
  }

  float TPhiMatrix::get_value_unsafe(const Token& token, int topic_id) const {
    auto iter = token2id_.find(token);

    if (iter == token2id_.end() || topic_id >= num_topics_) {
      throw std::runtime_error("TPhiMatrix::get_value: index outside bounds or no such token in matrix");
    }
    return data_[iter->second][topic_id];
  }

  void TPhiMatrix::set_value_unsafe(const Token& token, int topic_id, float value) {
    auto iter = token2id_.find(token);

    if (iter == token2id_.end() || topic_id >= num_topics_) {
      throw std::runtime_error("TPhiMatrix::set_value: index outside bounds or no such token in matrix");
    }
    data_[iter->second][topic_id] = value;
  }

}  // namespace ParallelEM
