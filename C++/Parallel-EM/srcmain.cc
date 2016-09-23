#include "common.h"
#include "matrix.h"
#include "processor.h"

using namespace ParallelEM;

namespace po = boost::program_options;
namespace bf = boost::filesystem;

namespace {
  struct TParameters {
    int num_topics;
    int num_outer_iter;
    int num_inner_iter;
    int num_processors;
    int num_top_tokens;
    std::string data_path;
  };

  void ParseAndPrintArgs(int argc, char* argv[], TParameters* p) {
    po::options_description all_options("PEM options");
    all_options.add_options()
      ("help,H", "Show help")
      ("num-topics,T", po::value(&p->num_topics)->default_value(DEFAULT_NUM_TOPICS), "Input number of topics")
      ("num-outer-iter,O", po::value(&p->num_outer_iter)->default_value(DEFAULT_NUM_O_ITER), "Input number of collection passes")
      ("num-inner-iter,I", po::value(&p->num_inner_iter)->default_value(DEFAULT_NUM_I_ITER), "Input number of document passes")
      ("num-top-tokens,K", po::value(&p->num_top_tokens)->default_value(DEFAULT_NUM_TOP_TOKENS), "Input number of top tokens")
      ("num-threads,P", po::value(&p->num_processors)->default_value(DEFAULT_NUM_PROCESSORS), "Input number of threads")
      ("data-path,D", po::value(&p->data_path)->default_value(DEFAULT_DATA_PATH), "Input path to files with documents")
      ;

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(all_options).run(), vm);
    notify(vm);

    std::cout << "num-topics:     " << p->num_topics << std::endl;
    std::cout << "num-outer-iter: " << p->num_outer_iter << std::endl;
    std::cout << "num-inner-iter: " << p->num_inner_iter << std::endl;
    std::cout << "num-top-tokens: " << p->num_top_tokens << std::endl;
    std::cout << "num-threads:    " << p->num_processors << std::endl;
    std::cout << "data-path:      " << p->data_path << std::endl;
  }

  void ReadFilenames(const std::string& data_path, std::shared_ptr<std::vector<std::string> > filenames) {
    bf::path data_path_var(data_path);
    bf::directory_iterator end_iter;

    if (!bf::exists(data_path_var) || !bf::is_directory(data_path_var)) {
      throw std::runtime_error("ReadFilenames: given data-path name is incorrect or not exists");
    }

    bool has_any_document = false;
    for (bf::directory_iterator dir_iter(data_path_var);
      dir_iter != end_iter; ++dir_iter) {
      if (bf::is_regular_file(*dir_iter) && dir_iter->path().extension() == ".txt") {
        has_any_document = true;
        filenames->push_back((data_path_var / dir_iter->path().filename()).string());
      }
    }
    if (!has_any_document) {
      throw std::runtime_error("ReadFilenames: given data-path doesn't contain any .txt file");
    }
  }

  void CreateDictionary(std::shared_ptr<std::vector<std::string> > filenames, const std::string& data_path,
                        std::vector<Token>* tokens, std::vector<float>* token_freqs) {
    std::unordered_map<Token, int> token_freq_map;
    bool loaded = false;
    // try to read dictionary from disk
    bf::path data_path_var(data_path);
    bf::directory_iterator end_iter;
    for (bf::directory_iterator dir_iter(data_path_var);
      dir_iter != end_iter; ++dir_iter) {
      if (bf::is_regular_file(*dir_iter) && dir_iter->path().extension() == ".dict") {
        std::cout << "\nStart loading dictionary..." << std::endl;
        loaded = true;
        std::ifstream file_stream;
        file_stream.open((data_path_var / dir_iter->path().filename()).string().c_str());

        while (!file_stream.eof()) {
          std::string str;
          std::getline(file_stream, str);

        if (str.empty())
          continue;

          std::vector<std::string> token_freq{ std::istream_iterator<std::string>{std::istringstream(str)},
            std::istream_iterator<std::string>{} };

          token_freq_map.insert(std::make_pair(token_freq[0], std::stoi(token_freq[1])));
        }
      }
    }

    if (!loaded) {
      std::cout << "\nStart gathering dictionary..." << std::endl;
      for (const auto& filename : *filenames) {
        std::ifstream file_stream;
        file_stream.open(filename.c_str());

        if (file_stream.fail()) {
          throw std::runtime_error("CreateDictionary: fail to open file " + filename);
        }

        while (!file_stream.eof()) {
          std::string str;
          std::getline(file_stream, str);

          if (str.empty())
            continue;

          std::vector<std::string> doc_tokens{ std::istream_iterator<std::string>{std::istringstream(str)},
            std::istream_iterator<std::string>{} };

          for (auto& t : doc_tokens) {
            auto iter = token_freq_map.find(t);
            if (iter == token_freq_map.end()) {
              token_freq_map.insert(std::make_pair(t, 1));
            } else {
              iter->second++;
            }
          }
        }
      }

      // save gathered dictionary to disk for future usage
      std::ofstream file_stream;
      file_stream.open((data_path_var / DICTIONARY_NAME).string());
      for (auto iter = token_freq_map.begin(); iter != token_freq_map.end(); ++iter) {
        file_stream << iter->first << " " << iter->second << "\n";
      }
      file_stream.close();
    }

    for (auto iter = token_freq_map.begin(); iter != token_freq_map.end(); ++iter) {
      tokens->push_back(iter->first);
      token_freqs->push_back(iter->second);
    }
    std::cout << "Complete creating dictionary\n" << std::endl;
  }

  void SaveTopTokens(std::shared_ptr<TPhiMatrix> phi, int num_top_tokens,
                     const std::vector<Token>& tokens, const std::string& data_path) {
    bf::path data_path_var(data_path);
    std::ofstream file_stream;
    file_stream.open((data_path_var / OUTPUT_FILENAME).string());
    for (int topic_id = 0; topic_id < phi->num_topics(); ++topic_id) {
      file_stream << "Topic_" << topic_id << " ";
      std::vector<std::pair<float, int> > topic_values;
      for (int token_id = 0; token_id < phi->num_tokens(); ++token_id) {
        topic_values.push_back(std::make_pair(phi->get_value_unsafe(tokens[token_id], topic_id), token_id));
      }

      std::sort(topic_values.begin(), topic_values.end());
      for (int token_id = 0; token_id < num_top_tokens; ++token_id) {
        file_stream << tokens[topic_values[topic_values.size() - token_id - 1].second] << " ("
                    << topic_values[topic_values.size() - token_id - 1].first << ") ";
      }
      file_stream << "\n";
    }
    file_stream.close();
  }

  bool HasBusyProcessor(const std::vector<std::shared_ptr<TProcessor> >& processors) {
    for (auto p : processors) {
      if (!p->is_stopped()) return true;
    }
    return false;
  }
}

int main(int argc, char* argv[]) {
  try {
    auto time_start = boost::posix_time::microsec_clock::local_time();
    TParameters p;

    ParseAndPrintArgs(argc, argv, &p);

    std::shared_ptr<std::vector<std::string> > filenames(new std::vector<std::string>);
    ReadFilenames(p.data_path, filenames);
    std::vector<Token> tokens;
    std::vector<float> token_freqs;
    CreateDictionary(filenames, p.data_path, &tokens, &token_freqs);
    auto time_finish = boost::posix_time::microsec_clock::local_time();
    std::cout << "\nElapsed time: " << (time_finish - time_start).total_seconds() << " sec.\n";

    std::shared_ptr<TPhiMatrix> phi(new TPhiMatrix(p.num_topics, tokens, true));

    float n = 0.0f;
    for (float value : token_freqs) {
      n += value;
    }

    time_start = boost::posix_time::microsec_clock::local_time();
    for (int outer_iter = 0; outer_iter < p.num_outer_iter; ++outer_iter) {
      std::shared_ptr<TPhiMatrix> n_wt(new TPhiMatrix(p.num_topics, tokens, false));
      std::shared_ptr<TPhiMatrix> n_t(new TPhiMatrix(p.num_topics, { TEMP_TOKEN }, false));
      std::shared_ptr<TFilenameQueue> filename_queue(new TFilenameQueue(*filenames));

      // run processors to proceed E-step and collect n_wt and n_t
      std::vector<float> perplexity_vector(p.num_processors, 0.0f);
      std::vector<std::shared_ptr<TProcessor> > processors;
      for (int j = 0; j < p.num_processors; ++j) {
        processors.push_back(std::shared_ptr<TProcessor>(
          new TProcessor(filename_queue, phi, n_wt, n_t, &(perplexity_vector[j]), p.num_inner_iter)));
      }

      while (HasBusyProcessor(processors)) {
        std::this_thread::sleep_for(std::chrono::seconds(CHECK_UPDATES_PERIOD));
      }

      // proceed M-step and update phi
      std::cout << std::endl;
      for (const auto& token : tokens) {
        for (int topic_id = 0; topic_id < p.num_topics; ++topic_id) {
          float denominator = n_t->get_value_unsafe(TEMP_TOKEN, topic_id);
          float value = denominator > 0.0f ? n_wt->get_value_unsafe(token, topic_id) / denominator : 0.0f;
          phi->set_value_unsafe(token, topic_id, value);
        }
      }

      // print perplexity and info
      float perplexity_value = 0.0f;
      for (float value : perplexity_vector) {
        perplexity_value += value;
      }

      time_finish = boost::posix_time::microsec_clock::local_time();
      std::cout << "Iter# " << outer_iter
                << ", Elapsed time: " << (time_finish - time_start).total_seconds()
                << ", Perplexity: " << exp(-(1.0f / n) * perplexity_value) << "\n";
    }
    time_finish = boost::posix_time::microsec_clock::local_time();
    std::cout << "\nElapsed time: " << (time_finish - time_start).total_seconds() << " sec.\n";
    std::cout << "\nTop tokens were saved into file.\n";
    SaveTopTokens(phi, p.num_top_tokens, tokens, p.data_path);
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
