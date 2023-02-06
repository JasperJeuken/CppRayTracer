/*
Bounding volume hierachy node. Used to optimize hit calculations.
*/
#ifndef BVH_H
#define BVH_H

#include <typeinfo>

#include "utility.h"
#include "hittable.h"
#include "hittable_list.h"

// Node in the BVH tree, contains a 'left' and 'right' object, and an axis-aligned bounding box
class bvh_node : public hittable {
public:
    bvh_node();

    bvh_node(const hittable_list& list, double time0 = 0.0, double time1 = 1.0) : bvh_node(list.objects, 0, list.objects.size(), time0, time1) {}

    bvh_node(const std::vector<shared_ptr<hittable>>& src_objects, size_t start, size_t end, double time0, double time1);

    // Override hit function
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

    // Override bounding box function
    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;

    // Override load normal map to share same image_texture for all objects in the tree
    virtual void load_normal_map(const char* filename, double strength = 1.0) override {
        normal_map = make_shared<image_texture>(filename);
        left->set_normal_map(normal_map);
        right->set_normal_map(normal_map);
    }

    // Override set normal map to share same image texture for all objects in the tree
    virtual void set_normal_map(shared_ptr<image_texture> map, double strength = 1.0) override {
        normal_map = map;
        has_normal_map = true;
        normal_map_strength = strength;
        left->set_normal_map(map, strength);
        right->set_normal_map(map, strength);
    }

public:
    shared_ptr<hittable> left;
    shared_ptr<hittable> right;
    aabb box;
};

// Compare two hittable objects along an axis
inline bool box_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b, int axis) {
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
        std::cout << "No bounding box in bvh_node constructor.\n";

    return box_a.min().e[axis] < box_b.min().e[axis];
}

// Compare on x-axis
bool box_x_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 0);
}

// Compare on y-axis
bool box_y_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 1);
}

// Compare on z-axis
bool box_z_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 2);
}

// Constructor that creates the BVH tree from a list of hittable objects
bvh_node::bvh_node(const std::vector<shared_ptr<hittable>>& src_objects, size_t start, size_t end, double time0, double time1
) {
    // Create a modifiable copy of the objects
    auto objects = src_objects;

    // Choose a random axis along which the BVH node is constructed
    int axis = random_int(0, 2);
    auto comparator = (axis == 0) ? box_x_compare
        : (axis == 1) ? box_y_compare
        : box_z_compare;

    size_t object_span = end - start;
    
    // If there is only one object, set both left and right to that object
    if (object_span == 1) {
        left = right = objects[start];
    }

    // If there are two objects, compare them and set to left and right
    else if (object_span == 2) {
        if (comparator(objects[start], objects[start + 1])) {
            left = objects[start];
            right = objects[start + 1];
        }
        else {
            left = objects[start + 1];
            right = objects[start];
        }
    }

    // If there are more than two objects, sort the objects and recursively make a BVH node for both left and right
    else {
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        auto mid = start + object_span / 2;
        left = make_shared<bvh_node>(objects, start, mid, time0, time1);
        right = make_shared<bvh_node>(objects, mid, end, time0, time1);
    }

    // Create a bounding box for this BVH node
    aabb box_left, box_right;
    if (!left->bounding_box(time0, time1, box_left)
        || !right->bounding_box(time0, time1, box_right)
        )
        std::cout << "No bounding box in bvh_node constructor.\n";
    box = surrounding_box(box_left, box_right);
}

// Determine if the BVH node is hit by a ray
bool bvh_node::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    if (box.hit(r, t_min, t_max)) {
        if (left->hit(r, t_min, t_max, rec)) {
            right->hit(r, t_min, rec.t, rec);
            return true;
        }
        else {
            return right->hit(r, t_min, t_max, rec);
        }
    }
    return false;
}

// Return the bounding box of the BVH node
bool bvh_node::bounding_box(double time0, double time1, aabb& output_box) const {
    output_box = box;
    return true;
}

#endif // !BVH_H
