#if 0
g++ -std=c++11 -shared concaveman.cpp -o libconcaveman.so
exit 0
#endif

//
// Author: Stanislaw Adaszewski, 2019
//

#include <memory>
#include <stdexcept>
#include <list>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <limits>
#include <set>
#include <queue>
#include <assert.h>


template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}


template<class T> class compare_first {
public:
    bool operator()(const T &a, const T &b) {
        return (std::get<0>(a) < std::get<0>(b));
    }
};


template<class T> T orient2d(
    const std::array<T, 2> &p1,
    const std::array<T, 2> &p2,
    const std::array<T, 2> &p3) {

    T res = (p2[1] - p1[1]) * (p3[0] - p2[0]) -
        (p2[0] - p1[0]) * (p3[1] - p2[1]);

    return res;
}


template<class T> bool intersects(
    const std::array<T, 2> &p1,
    const std::array<T, 2> &q1,
    const std::array<T, 2> &p2,
    const std::array<T, 2> &q2) {

    auto res = (p1[0] != q2[0] || p1[1] != q2[1]) &&
        (q1[0] != p2[0] || q1[1] != p2[1]) &&
        (orient2d(p1, q1, p2) > 0) != (orient2d(p1, q1, q2) > 0) &&
        (orient2d(p2, q2, p1) > 0) != (orient2d(p2, q2, q1) > 0);

    return res;
}


template<class T> T getSqDist(
    const std::array<T, 2> &p1,
    const std::array<T, 2> &p2) {

    auto dx = p1[0] - p2[0];
    auto dy = p1[1] - p2[1];
    return dx * dx + dy * dy;
}


template<class T> T sqSegDist(
    const std::array<T, 2> &p,
    const std::array<T, 2> &p1,
    const std::array<T, 2> &p2) {

    auto x = p1[0];
    auto y = p1[1];
    auto dx = p2[0] - x;
    auto dy = p2[1] - y;

    if (dx != 0 || dy != 0) {
        auto t = ((p[0] - x) * dx + (p[1] - y) * dy) / (dx * dx + dy * dy);
        if (t > 1) {
            x = p2[0];
            y = p2[1];
        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
        }
    }

    dx = p[0] - x;
    dy = p[1] - y;

    return dx * dx + dy * dy;
}


template<class T> T sqSegSegDist(T x0, T y0,
    T x1, T y1,
    T x2, T y2,
    T x3, T y3) {
    auto ux = x1 - x0;
    auto uy = y1 - y0;
    auto vx = x3 - x2;
    auto vy = y3 - y2;
    auto wx = x0 - x2;
    auto wy = y0 - y2;
    auto a = ux * ux + uy * uy;
    auto b = ux * vx + uy * vy;
    auto c = vx * vx + vy * vy;
    auto d = ux * wx + uy * wy;
    auto e = vx * wx + vy * wy;
    auto D = a * c - b * b;

    T sc, sN, tc, tN;
    auto sD = D;
    auto tD = D;

    if (D == 0) {
        sN = 0;
        sD = 1;
        tN = e;
        tD = c;
    } else {
        sN = b * e - c * d;
        tN = a * e - b * d;
        if (sN < 0) {
            sN = 0;
            tN = e;
            tD = c;
        } else if (sN > sD) {
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0) {
        tN = 0;
        if (-d < 0) sN = 0;
        else if (-d > a) sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    } else if (tN > tD) {
        tN = tD;
        if (-d + b < 0) sN = 0;
        else if (-d + b > a) sN = sD;
        else {
            sN = -d + b;
            sD = a;
        }
    }

    sc = ((sN == 0) ? 0 : sN / sD);
    tc = ((tN == 0) ? 0 : tN / tD);

    auto cx = (1 - sc) * x0 + sc * x1;
    auto cy = (1 - sc) * y0 + sc * y1;
    auto cx2 = (1 - tc) * x2 + tc * x3;
    auto cy2 = (1 - tc) * y2 + tc * y3;
    auto dx = cx2 - cx;
    auto dy = cy2 - cy;

    return dx * dx + dy * dy;
}


template<class T, int DIM, int MAX_CHILDREN, class DATA> class rtree {
public:
    typedef rtree<T, DIM, MAX_CHILDREN, DATA> type;
    typedef const type const_type;
    typedef type *type_ptr;
    typedef const type *type_const_ptr;
    typedef std::array<T, DIM * 2> bounds_type;
    typedef DATA data_type;

    rtree():
    m_data(), m_is_leaf(false) {
        for (auto i = 0; i < DIM; i++) {
            m_bounds[i] = std::numeric_limits<T>::max();
            m_bounds[i + DIM] = std::numeric_limits<T>::min();
        }
    }

    rtree(data_type data, const bounds_type &bounds):
    m_data(data), m_bounds(bounds), m_is_leaf(true) {
        for (auto i = 0; i < DIM; i++)
            if (bounds[i] > bounds[i + DIM])
                throw std::runtime_error("Bounds minima have to be less than maxima");
    }

    void insert(data_type data, const bounds_type &bounds) {
        if (m_is_leaf)
          throw std::runtime_error("Cannot insert into leaves");

        m_bounds = updated_bounds(bounds);
        if (m_children.size() < MAX_CHILDREN) {
            auto r = make_unique<type>(data, bounds);
            m_children.push_back(std::move(r));
            return;
        }

        std::reference_wrapper<type> best_child = *m_children.begin()->get();
        auto best_volume = volume(best_child.get().updated_bounds(bounds));
        for (auto it = ++m_children.begin(); it != m_children.end(); it++) {
            auto v = volume((*it)->updated_bounds(bounds));
            if (v < best_volume) {
                best_volume = v;
                best_child = *it->get();
            }
        }
        // std::cout << "best_child: " << best_child.get().id() << std::endl;
        if (!best_child.get().is_leaf()) {
            best_child.get().insert(data, bounds);
            return;
        }

        auto leaf = make_unique<type>(best_child.get().data(),
            best_child.get().bounds());
        best_child.get().m_is_leaf = false;
        best_child.get().m_data = data_type();
        best_child.get().m_children.push_back(std::move(leaf));
        best_child.get().insert(data, bounds);
    }

    void intersection(const bounds_type &bounds,
        std::vector<std::reference_wrapper<const_type>> &res) const {
        if (!intersects(bounds))
            return;
        if (m_is_leaf) {
            res.push_back(*this);
            return;
        }
        for (auto &ch : m_children)
            ch->intersection(bounds, res);
    }

    std::vector<std::reference_wrapper<const_type>> intersection(const bounds_type& bounds) const {
        std::vector<std::reference_wrapper<const_type>> res;
        intersection(bounds, res);
        return res;
    }

    bool intersects(const bounds_type &bounds) const {
        for (auto i = 0; i < DIM; i++) {
            if (m_bounds[i] > bounds[i + DIM])
                return false;
            if (m_bounds[i + DIM] < bounds[i])
                return false;
        }
        return true;
    }

    void erase(data_type data, const bounds_type &bounds) {
        if (m_is_leaf)
            throw std::runtime_error("Cannot erase from leaves");

        if (!intersects(bounds))
            return;

        for (auto it = m_children.begin(); it != m_children.end(); ) {
            if (!(*it)->m_is_leaf) {
                (*it)->erase(data, bounds);
                it++;
            } else if ((*it)->m_data == data &&
                (*it)->m_bounds == bounds) {
                m_children.erase(it++);
            } else
                it++;
        }
    }

    bounds_type updated_bounds(const bounds_type &child_bounds) const {
        bounds_type res;
        for (auto i = 0; i < DIM; i++) {
            res[i] = std::min(child_bounds[i], m_bounds[i]);
            res[i + DIM] = std::max(child_bounds[i + DIM], m_bounds[i + DIM]);
        }
        return res;
    }

    static T volume(const bounds_type &bounds) {
        T res = 1;
        for (auto i = 0; i < DIM; i++) {
            auto delta = bounds[i + DIM] - bounds[i];
            res *= delta;
        }
        return res;
    }

    const bounds_type& bounds() const {
        return m_bounds;
    }

    bool is_leaf() const {
        return m_is_leaf;
    }

    data_type data() const {
        return m_data;
    }

    const std::list<std::unique_ptr<type>>& children() const {
        return m_children;
    }

    static std::string bounds_to_string(const bounds_type &bounds) {
        std::string res = "( ";
        for (auto i = 0; i < DIM * 2; i++) {
            if (i > 0)
                res += ", ";
            res += std::to_string(bounds[i]);
        }
        res += " )";
        return res;
    }

    void to_string(std::string &res, int tab) const {
        std::string pad(tab, '\t');

        if (m_is_leaf) {
            res += pad + "{ data: " + std::to_string(m_data) +
                ", bounds: " + bounds_to_string(m_bounds) +
                " }";
            return;
        }

        res += pad + "{ bounds: " + bounds_to_string(m_bounds) +
            ", children: [\n";
        auto i = 0;
        for (auto &ch : m_children) {
            if (i++ > 0)
                res += "\n";
            ch->to_string(res, tab + 1);
        }
        res += "\n" + pad + "]}";
    }

    std::string to_string() const {
        std::string res;
        to_string(res, 0);
        return res;
    }

private:
    bool m_is_leaf;
    data_type m_data;
    std::list<std::unique_ptr<type>> m_children;
    bounds_type m_bounds;
};


template<class T> struct Node {
    typedef Node<T> type;
    typedef type *type_ptr;
    typedef std::array<T, 2> point_type;

    Node(): p(),
    minX(), minY(), maxX(), maxY() {

    }

    Node(const point_type &p): Node() {
        this->p = p;
    }

    point_type p;
    T minX;
    T minY;
    T maxX;
    T maxY;
};


template <class T> class CircularList;


template<class T> class CircularElement {
public:
    typedef CircularElement<T> type;
    typedef type *ptr_type;

    template<class... Args> CircularElement<T>(Args&&... args):
    m_data(std::forward<Args>(args)...) {

    }

    T& data() {
        return m_data;
    }

    template<class... Args> CircularElement<T>* insert(Args&&... args) {
        auto elem = new CircularElement<T>(std::forward<Args>(args)...);
        elem->m_prev = this;
        elem->m_next = m_next;
        m_next->m_prev = elem;
        m_next = elem;
        return elem;
    }

    CircularElement<T>* prev() {
        return m_prev;
    }

    CircularElement<T>* next() {
        return m_next;
    }

private:
    T m_data;
    CircularElement<T> *m_prev;
    CircularElement<T> *m_next;

    friend class CircularList<T>;
};


template<class T> class CircularList {
public:
    typedef CircularElement<T> element_type;

    CircularList(): m_last(nullptr) {

    }

    ~CircularList() {
        std::cout << "~CircularList()" << std::endl;
        auto node = m_last;
        auto i = 0;
        while (true) {
            // std::cout << (i++) << std::endl;
            auto tmp = node;
            node = node->m_next;
            delete tmp;
            if (node == m_last)
                break;
        }
    }

    template<class... Args> CircularElement<T>* insert(element_type *prev, Args&&... args) {
        auto elem = new CircularElement<T>(std::forward<Args>(args)...);

        if (prev == nullptr && m_last != nullptr)
            throw std::runtime_error("Once the list is non-empty you must specify where to insert");

        if (prev == nullptr) {
            elem->m_prev = elem->m_next = elem;
        } else {
            elem->m_prev = prev;
            elem->m_next = prev->m_next;
            prev->m_next->m_prev = elem;
            prev->m_next = elem;
        }

        m_last = elem;

        return elem;
    }


private:
    element_type *m_last;
};


template<class T> void updateBBox(typename CircularElement<T>::ptr_type elem) {
    auto &node(elem->data());
    auto p1 = node.p;
    auto p2 = elem->next()->data().p;
    node.minX = std::min(p1[0], p2[0]);
    node.minY = std::min(p1[1], p2[1]);
    node.maxX = std::max(p1[0], p2[0]);
    node.maxY = std::max(p1[1], p2[1]);
}


template<class T, int MAX_CHILDREN> std::vector<std::array<T, 2>> concaveman(
    const std::vector<std::array<T, 2>> &points,
    const std::vector<int> &hull,
    T concavity=2, T lengthThreshold=0) {

    typedef Node<T> node_type;
    typedef std::array<T, 2> point_type;
    typedef CircularElement<node_type> circ_elem_type;
    typedef CircularList<node_type> circ_list_type;
    typedef circ_elem_type *circ_elem_ptr_type;

    std::cout << "concaveman()" << std::endl;

    if (hull.size() == points.size()) {
        std::vector<point_type> res;
        for (auto &i : hull) res.push_back(points[i]);
        return res;
    }

    rtree<T, 2, MAX_CHILDREN, point_type> tree;

    for (auto &p : points)
        tree.insert(p, { p[0], p[1], p[0], p[1] });

    circ_list_type circList;
    circ_elem_ptr_type last = nullptr;

    std::list<circ_elem_ptr_type> queue;

    for (auto &idx : hull) {
        auto &p = points[idx];
        tree.erase(p, { p[0], p[1], p[0], p[1] });
        last = circList.insert(last, p);
        queue.push_back(last);
    }

    std::cout << "Starting hull: ";
    for (auto elem = last->next(); ; elem=elem->next()) {
        std::cout << elem->data().p[0] << " " << elem->data().p[1] << std::endl;
        if (elem == last)
            break;
    }

    rtree<T, 2, MAX_CHILDREN, circ_elem_ptr_type> segTree;
    for (auto &elem : queue) {
        auto &node(elem->data());
        updateBBox<node_type>(elem);
        segTree.insert(elem, { node.minX,
            node.minY, node.maxX, node.maxY });
    }

    auto sqConcavity = concavity * concavity;
    auto sqLenThreshold = lengthThreshold * lengthThreshold;

    while (!queue.empty()) {
        auto elem = *queue.begin();
        queue.pop_front();

        auto a = elem->data().p;
        auto b = elem->next()->data().p;

        auto sqLen = getSqDist(a, b);
        if (sqLen < sqLenThreshold)
            continue;

        auto maxSqLen = sqLen / sqConcavity;

        bool ok;
        auto p = findCandidate(tree, elem->prev()->data().p, a, b,
            elem->next()->next()->data().p, maxSqLen, segTree, ok);

        if (ok && std::min(getSqDist(p, a), getSqDist(p, b)) <= maxSqLen) {

            std::cout << "Modifying hull, p: " << p[0] << " " << p[1] << std::endl;

            queue.push_back(elem);
            queue.push_back(elem->insert(p));

            auto &node = elem->data();
            auto &next = elem->next()->data();

            tree.erase(p, { p[0], p[1], p[0], p[1] });
            segTree.erase(elem, { node.minX, node.minY, node.maxX, node.maxY });

            updateBBox<node_type>(elem);
            updateBBox<node_type>(elem->next());

            segTree.insert(elem, { node.minX, node.minY, node.maxX, node.maxY });
            segTree.insert(elem->next(), { next.minX, next.minY, next.maxX, next.maxY });
        }
    }

    std::vector<point_type> concave;
    for (auto elem = last->next(); ; elem = elem->next()) {
        concave.push_back(elem->data().p);
        if (elem == last)
            break;
    }

    return concave;
}


template<class T, int MAX_CHILDREN> std::array<T, 2> findCandidate(
    const rtree<T, 2, MAX_CHILDREN, std::array<T, 2>> &tree,
    const std::array<T, 2> &a,
    const std::array<T, 2> &b,
    const std::array<T, 2> &c,
    const std::array<T, 2> &d,
    T maxDist,
    const rtree<T, 2, MAX_CHILDREN, typename CircularElement<Node<T>>::ptr_type> &segTree,
    bool &ok) {

    typedef std::array<T, 2> point_type;
    typedef CircularElement<Node<T>> circ_elem_type;
    typedef rtree<T, 2, MAX_CHILDREN, std::array<T, 2>> tree_type;
    typedef const tree_type const_tree_type;
    typedef std::reference_wrapper<const_tree_type> tree_ref_type;
    typedef std::tuple<T, tree_ref_type> tuple_type;

    std::cout << "findCandidate(), maxDist: " << maxDist << std::endl;

    ok = false;

    std::priority_queue<tuple_type, std::vector<tuple_type>, compare_first<tuple_type>> queue;
    std::reference_wrapper<const_tree_type> node = tree;

    while (true) {
        for (auto &child : node.get().children()) {

            auto bounds = child->bounds();
            point_type pt = { bounds[0], bounds[1] };

            auto dist = child->is_leaf() ? sqSegDist(pt, b, c) : sqSegBoxDist(b, c, *child);
            if (dist > maxDist)
                continue;

            queue.push(tuple_type(-dist, *child));
        }

        while (!queue.empty() && std::get<1>(queue.top()).get().is_leaf()) {
            auto item = queue.top();
            queue.pop();

            auto bounds = std::get<1>(item).get().bounds();
            point_type p = { bounds[0], bounds[1] };

            auto d0 = sqSegDist(p, a, b);
            auto d1 = sqSegDist(p, c, d);

            if (std::get<0>(item) < d0 && std::get<0>(item) < d1 &&
                noIntersections(b, p, segTree) &&
                noIntersections(c, p, segTree)) {

                ok = true;
                return std::get<1>(item).get().data();
            }
        }

        if (queue.empty())
            break;

        node = std::get<1>(queue.top());
        queue.pop();
    }

    return point_type();
}


template<class T, int MAX_CHILDREN, class USER_DATA> T sqSegBoxDist(
    const std::array<T, 2> &a,
    const std::array<T, 2> &b,
    const rtree<T, 2, MAX_CHILDREN, USER_DATA> &bbox) {

    if (inside(a, bbox) || inside(b, bbox))
        return 0;

    auto &bounds = bbox.bounds();
    auto minX = bounds[0];
    auto minY = bounds[1];
    auto maxX = bounds[2];
    auto maxY = bounds[3];

    auto d1 = sqSegSegDist(a[0], a[1], b[0], b[1], minX, minY, maxX, minY);
    if (d1 == 0) return 0;

    auto d2 = sqSegSegDist(a[0], a[1], b[0], b[1], minX, minY, minX, maxY);
    if (d2 == 0) return 0;

    auto d3 = sqSegSegDist(a[0], a[1], b[0], b[1], maxX, minY, maxX, maxY);
    if (d3 == 0) return 0;

    auto d4 = sqSegSegDist(a[0], a[1], b[0], b[1], minX, maxY, maxX, maxY);
    if (d4 == 0) return 0;

    return std::min(std::min(d1, d2), std::min(d3, d4));
}


template<class T, int MAX_CHILDREN, class USER_DATA> bool inside(
    const std::array<T, 2> &a,
    const rtree<T, 2, MAX_CHILDREN, USER_DATA> &bbox) {

    auto &bounds = bbox.bounds();

    auto minX = bounds[0];
    auto minY = bounds[1];
    auto maxX = bounds[2];
    auto maxY = bounds[3];

    auto res = (a[0] >= minX) &&
        (a[0] <= maxX) &&
        (a[1] >= minY) &&
        (a[1] <= maxY);
    return res;
}


template<class T, int MAX_CHILDREN> bool noIntersections(
    const std::array<T, 2> &a,
    const std::array<T, 2> &b,
    const rtree<T, 2, MAX_CHILDREN, typename CircularElement<Node<T>>::ptr_type> &segTree) {

    auto minX = std::min(a[0], b[0]);
    auto minY = std::min(a[1], b[1]);
    auto maxX = std::max(a[0], b[0]);
    auto maxY = std::max(a[1], b[1]);

    auto isect = segTree.intersection({ minX, minY, maxX, maxY });

    for (decltype(segTree) &ch : isect) {
        auto elem = ch.data();

        if (intersects(elem->data().p, elem->next()->data().p, a, b))
            return false;
    }

    return true;
}


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

extern "C" {
    void pyconcaveman2d(double *points_c, size_t num_points,
        int *hull_points_c, size_t num_hull_points,
        double concavity, double lengthThreshold,
        double **concave_points_c, size_t *num_concave_points,
        void(**p_free)(void*));
}

void pyconcaveman2d(double *points_c, size_t num_points,
    int *hull_points_c, size_t num_hull_points,
    double concavity, double lengthThreshold,
    double **p_concave_points_c,
    size_t *p_num_concave_points,
    void(**p_free)(void*)) {

    std::cout << "pyconcaveman2d(), concavity: " << concavity <<
        " lengthThreshold: " << lengthThreshold << std::endl;

    typedef double T;
    typedef std::array<T, 2> point_type;

    std::vector<point_type> points(num_points);
    for (auto i = 0; i < num_points; i++) {
        points[i] = { points_c[i << 1], points_c[(i << 1) + 1] };
    }

    std::cout << "points:" << std::endl;
    for (auto &p : points)
        std::cout << p[0] << " " << p[1] << std::endl;

    std::vector<int> hull(num_hull_points);
    for (auto i = 0; i < num_hull_points; i++) {
        hull[i] = hull_points_c[i];
    }

    std::cout << "hull:" << std::endl;
    for (auto &i : hull)
        std::cout << i << std::endl;

    auto concave_points = concaveman<T, 16>(points, hull, concavity, lengthThreshold);

    std::cout << "concave_points:" << std::endl;
    for (auto &p : concave_points)
        std::cout << p[0] << " " << p[1] << std::endl;

    double *concave_points_c = *p_concave_points_c = (double*) malloc(sizeof(double) * 2 * concave_points.size());
    for (auto i = 0; i < concave_points.size(); i++) {
        concave_points_c[i << 1] = concave_points[i][0];
        concave_points_c[(i << 1) + 1] = concave_points[i][1];
    }

    *p_num_concave_points = concave_points.size();
    *p_free = free;
}
