#if 0
g++ -std=c++11 -I../../../main/cpp -o test_concaveman
exit 0
#endif

static void test_00_rtree_basic() {
    typedef rtree<float, 2, 8, intptr_t> myrtree;
    typedef rtree<float, 2, 16, intptr_t> myrtree2;
    typedef rtree<float, 2, 5, intptr_t> myrtree3;
    typedef rtree<float, 2, 32, intptr_t> myrtree4;
    myrtree tree;
    myrtree2 tree2;
    myrtree3 tree3;
    myrtree4 tree4;

    // myrtree::bounds_type bounds ;
    for (auto i = 0; i < 500; i++) {
        auto a = static_cast<float>(random() % 100);
        auto b = static_cast<float>(random() % 100);
        auto c = static_cast<float>(random() % 100);
        auto d = static_cast<float>(random() % 100);
        myrtree::bounds_type bounds = { std::min(a, c), std::min(b, d),
            std::max(a, c), std::max(b, d) };
        tree.insert(i, bounds);
        tree2.insert(i, bounds);
        tree3.insert(i, bounds);
        tree4.insert(i, bounds);
        // std::cout << "i: " << i << ", tree: " << tree.to_string() << std::endl;
    }
    // tree.insert(1, { 1, 2, 1, 2 });
    // tree.insert(2, { 0, 1, 0, 1 });
    // tree.insert(2, { 0, 1, 0, 3 });
    auto isect = tree.intersection({0, 0, 50, 50});
    auto isect2 = tree2.intersection({0, 0, 50, 50});
    std::set<intptr_t> s, s2;
    for (auto &node : isect)
        s.insert(node.get().data());
    for (auto &node : isect2)
        s2.insert(node.get().data());
    std::cout << s.size() << " " << (s == s2) << std::endl;
    for (auto data : s)
        std::cout << data << " ";
    std::cout << std::endl;
    // std::cout << tree.to_string() << std::endl;
}

/* void test_01_rtree_erase() {
    typedef rtree<float, 2, 8> myrtree;
    myrtree tree;
    tree.insert(1, { 1, 2, 1, 2 });
    tree.insert(2, { 0, 1, 0, 1 });
    tree.insert(3, { 0, 1, 0, 3 });
    std::cout << tree.to_string() << std::endl << std::endl;
    tree.erase(1, { 1, 2, 1, 2 });
    tree.erase(2, { 0, 1, 0, 1 });
    std::cout << tree.to_string() << std::endl << std::endl;
    tree.insert(1, { 1, 2, 1, 2 });
    tree.insert(2, { 0, 1, 0, 1 });
    std::cout << tree.to_string() << std::endl << std::endl;
} */

/* void test_02_node_deleter() {
    std::cout << "test_02_node() : ";
    typedef double T;
    typedef Node<T> node_type;
    typename node_type::type_ptr last = nullptr;
    last = new node_type({ 5.0, 6.0 });
    last->insert({ 4.0, 3.0 });
    last->insert({ 1.0, 2.0 });
    last->free();
    delete last;
    std::cout << "PASSED" << std::endl;
} */

static void test_03_sqSegDist() {
    std::cout << "test_03_sqSegDist() : ";

    auto a = sqSegDist<double>({ 0, 0 }, { 0, 1 }, { 1, 0 });
    std::cout << "a: " << a << std::endl;
    assert(a == 0.5);

    auto b = sqSegDist<double>({ 0, 1 }, { 0, 1 }, { 1, 0 });
    std::cout << "b: " << b << std::endl;
    assert(b == 0);

    auto c = sqSegDist<double>({ -1, 0 }, { 0, 0 }, { 0, 1 });
    std::cout << "c: " << c << std::endl;
    assert(c == 1);

    auto d = sqSegDist<double>({ -1, -1 }, { 0, 0 }, { 0, 1 });
    std::cout << "d: " << d << std::endl;
    assert(d == 2);

    std::cout << "PASSED" << std::endl;
}


static void test_04_sqSegSegDist() {

}


static void test_05_CircularList() {
    std::cout << "test_05_CircularList() : ";
    typedef double T;
    typedef Node<T> node_type;
    typedef typename node_type::point_type point_type;
    auto lst = make_unique<CircularList<node_type>>();
    auto a = lst->insert(nullptr, point_type { 1, 2 });
    auto b = a->insert(point_type { 3, 4 });
    auto c = b->insert(point_type { 5, 6 });
    std::cout << a->data().p[0] << " " <<
        b->data().p[0] << " " << c->data().p[0] << std::endl;
    {
        std::unique_ptr<CircularList<node_type>> taker(std::move(lst));
    }
    assert(lst.get() == nullptr);
    std::cout << "PASSED" << std::endl;
}


static void test_06_priority_queue() {
    std::cout << "test_06_priority_queue() : ";
    typedef double T;
    typedef Node<T> node_type;
    typedef std::tuple<T, node_type> tuple_type;
    std::priority_queue<tuple_type, std::vector<tuple_type>, compare_first<tuple_type>> queue;

    queue.push(std::make_tuple(-5.0, node_type({ 1, 2 })));
    queue.push(std::make_tuple(-4.0, node_type({ 3, 4 })));
    queue.push(std::make_tuple(-1.0, node_type({ 5, 6 })));
    queue.push(std::make_tuple(-10.0, node_type({ 7, 8 })));
    queue.push(std::make_tuple(-0.5, node_type({ 9, 10 })));

    while (!queue.empty()) {
        std::cout << std::get<0>(queue.top()) << std::endl;
        queue.pop();
    }

    std::cout << "PASSED" << std::endl;
}


static void test_07_concaveman() {
    std::cout << "test_07_concaveman() : ";
    typedef double T;
    typedef std::array<T, 2> point_type;
    std::vector<point_type> points {
        {0, 0},
        {1, 0},
        {0.25, 0.15},
        {1, 1}
    };
    std::vector<int> hull {
        0, 1, 3
    };
    auto concave = concaveman<T, 16>(points, hull, 2, 1);
    for (auto &p : concave) {
        std::cout << p[0] << " " << p[1] << std::endl;
    }
    std::cout << "PASSED" << std::endl;
}


static void test_08_intersects() {
    std::cout << "test_08_intersects() : ";

    std::cout << intersects<double>({ 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 }) << std::endl;
    std::cout << intersects<double>({ 0, 0 }, { 1, 0 }, { 0, 0 }, { 1, 0 }) << std::endl;
    std::cout << intersects<double>({ 0, 0 }, { 1, 0 }, { 0.5, -0.5 }, { 0.5, 0.5 }) << std::endl;

    std::cout << "PASSED" << std::endl;
}


static void test_09_orient2d() {
    std::cout << "test_09_orient2d() : ";

    std::cout << orient2d<double>({ 0, 0 }, { 1, 0 }, { 0, 1 }) << std::endl;
    std::cout << orient2d<double>({ 0, 0 }, { 1, 0 }, { 0, 0 }) << std::endl;
    std::cout << orient2d<double>({ 0, 0 }, { 1, 0 }, { 0.5, -0.5 }) << std::endl;
    std::cout << orient2d<double>({ 0, 1 }, { 1, 0 }, { 0, 0 }) << std::endl;

    std::cout << "PASSED" << std::endl;
}

#if 0

int main() {
    // test_00_rtree_basic();
    // test_01_rtree_erase();
    // test_02_node_deleter();
    // test_03_sqSegDist();
    // test_05_CircularList();
    // test_06_priority_queue();
    test_07_concaveman();
    // test_08_intersects();
    // test_09_orient2d();
}

#endif
