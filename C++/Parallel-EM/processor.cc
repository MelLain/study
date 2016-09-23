#include "processor.h"

TProcessor::TProcessor(std::shared_ptr<TFilenameQueue> filename_queue,
  std::shared_ptr<TPhiMatrix> phi,
  std::shared_ptr<TPhiMatrix> n_wt,
  std::shared_ptr<TPhiMatrix> n_t,
  float* perplexity_value,
  int num_inner_iter)
    : filename_queue_(filename_queue)
    , phi_(phi)
    , n_wt_(n_wt)
    , n_t_(n_t)
    , perplexity_value_(perplexity_value)
    , num_inner_iter_(num_inner_iter)
    , is_stopped_(false) {
  boost::thread t(&TProcessor::ThreadFunction, this);
  thread_.swap(t);
}

TProcessor::~TProcessor() {
  is_stopped_ = true;
  if (thread_.joinable()) {
    thread_.join();
  }
}

void TProcessor::ThreadFunction() {
  std::string thread_id = boost::lexical_cast<std::string>(boost::this_thread::get_id());
  boost::mutex local_mutex;
  try {
    std::string current_filename;
    while (true) {
      {
        boost::lock_guard<boost::mutex> guard(local_mutex);
        if (filename_queue_->is_empty()) {
          break;
        }
        current_filename = filename_queue_->pop_back();
      }

      std::ifstream file_stream;
      file_stream.open((current_filename).c_str());

      if (file_stream.fail()) {
        boost::lock_guard<boost::mutex> guard(local_mutex);
        throw std::runtime_error("ThreadFunction: fail to open file " + current_filename);
      }

      while (!file_stream.eof()) {
        std::string str;
        std::getline(file_stream, str);

        if (str.empty())
          continue;

        std::vector<std::string> doc_tokens{ std::istream_iterator<std::string>{std::istringstream(str)},
          std::istream_iterator<std::string>{} };

        float n_d = 0;
        std::unordered_map<Token, int> n_dw;
        for (auto& t : doc_tokens) {
          auto iter = n_dw.find(t);
          if (iter == n_dw.end()) {
            n_dw.insert(std::make_pair(t, 1));
            ++n_d;
          }
          else {
            ++iter->second;
            ++n_d;
          }
        }

        std::vector<float> theta_d(phi_->num_topics(), 1.0f / phi_->num_topics());
        std::vector<float> Z_w;
        for (int inner_iter = 0; inner_iter < num_inner_iter_; ++inner_iter) {
          Z_w.clear();

          // count normalization constants
          for (auto iter = n_dw.begin(); iter != n_dw.end(); ++iter) {
            float z = 0.0f;
            for (int topic_id = 0; topic_id < phi_->num_topics(); ++topic_id) {
              z += phi_->get_value(iter->first, topic_id) * theta_d[topic_id];
            }
            Z_w.push_back(z);
          }

          // proceed E-step and update theta
          for (int topic_id = 0; topic_id < phi_->num_topics(); ++topic_id) {
            float value = 0.0f;
            int token_id = 0;
            for (auto iter = n_dw.begin(); iter != n_dw.end(); ++iter) {
              value += Z_w[token_id] > 0.0f ? (iter->second * phi_->get_value(iter->first, topic_id)
                * theta_d[topic_id]) / Z_w[token_id] : 0.0f;
              ++token_id;
            }
            theta_d[topic_id] = (1.0f / n_d) * value;
          }
        }

        // count perplexity
        for (auto iter = n_dw.begin(); iter != n_dw.end(); ++iter) {
          float value = 0.0f;
          for (int topic_id = 0; topic_id < phi_->num_topics(); ++topic_id) {
            value += phi_->get_value(iter->first, topic_id) * theta_d[topic_id];
          }
          *perplexity_value_ += iter->second * log(value > 0.0f ? value : 1.0f);
        }

        // update counters
        for (int topic_id = 0; topic_id < phi_->num_topics(); ++topic_id) {
          int token_id = 0;
          for (auto iter = n_dw.begin(); iter != n_dw.end(); ++iter) {
            float value = Z_w[token_id] > 0.0f ? (iter->second * phi_->get_value(iter->first, topic_id)
              * theta_d[topic_id]) / Z_w[token_id] : 0.0f;
            ++token_id;
            n_wt_->add_value(iter->first, topic_id, value);
            n_t_->add_value(TEMP_TOKEN, topic_id, value);
          }
        }
      }
    }
    is_stopped_ = true;
  }
  catch (std::exception& e) {
    boost::lock_guard<boost::mutex> guard(local_mutex);
    std::cout << thread_id << " - " << e.what() << std::endl;
    is_stopped_ = true;
  }
}
