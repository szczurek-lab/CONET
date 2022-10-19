
#include "../../src/tree/vertex_label_sampler.h"
#include "../test_utils.h"
std::vector<size_t> chromosome_markers{10, 20, 45};
size_t max_locus = 44;


void basic_operations_test() {
    BEGIN_TEST;
    VertexLabelSampler<double> sampler{max_locus, chromosome_markers};
    
    IS_TRUE(sampler.has_free_labels());

    sampler.add_label(std::make_pair(0, 1));
    sampler.remove_label(std::make_pair(0, 1));

    auto log_kernel = sampler.get_sample_label_log_kernel(); 
    auto expected_log_kernel = -std::log(45.0 + 45.0 + 300.0);
    IS_TRUE(std::abs(log_kernel - expected_log_kernel) <= 0.01);

    sampler.add_label(std::make_pair(0, 10));

    IS_FALSE(sampler.can_swap_one_breakpoint(std::make_pair(0, 5), std::make_pair(10, 11), 0, 0));
    END_TEST;
}


int main(void) {
    basic_operations_test();

}