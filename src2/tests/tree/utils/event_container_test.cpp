#include <iostream>
#include <set>
#include <random>

#include "../../src/utils/event_container.h"
#include "../test_utils.h"

using namespace::std;

void basic_operations_test()
{   
    BEGIN_TEST;

    EventContainer container = EventContainer(10);

    IS_TRUE(container.empty());
    
    Event ev1 = std::make_pair(0, 9);
    
    IS_FALSE(container.find(ev1));
    
    container.insert(ev1);
    
    IS_EQUAL(container.get_nth(0), ev1);
    
    Event ev2 = std::make_pair(0, 3);
    container.insert(ev2);
    
    IS_EQUAL(container.get_nth(0), ev2);
    IS_EQUAL(container.size(), 2);

    container.erase(ev2);
    IS_EQUAL(container.size(), 1);
    IS_EQUAL(container.get_nth(0), ev1);
    
    END_TEST;
}

Event sample_event(std::uniform_int_distribution<size_t> &dis, std::mt19937 &gen) {
    auto left = dis(gen);
    auto right = dis(gen); 
    if (left < right) {
        return std::make_pair(left, right);
    } else if (left > right) {
        return std::make_pair(right, left);
    } else {
        return sample_event(dis, gen);
    }
}

void randomized_test() {
    BEGIN_TEST;

    const size_t MAX_LOCUS = 10;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dis(0, MAX_LOCUS);

    EventContainer container = EventContainer(MAX_LOCUS);
    const size_t SAMPLES = 10000;
    std::set<Event> events; 

    for (size_t i = 0; i < SAMPLES; i++) {
        auto ev = sample_event(dis, gen); 
        if (events.find(ev) == events.end()) {
            events.insert(ev);
            container.insert(ev);
            IS_EQUAL(container.size(), events.size());

        } else {
            IS_TRUE(container.find(ev));
        }

    }

    END_TEST;
}

int main(void) {
    basic_operations_test();
    randomized_test();
}