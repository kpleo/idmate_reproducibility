"""Orthographic, vector ball-and-stick views of declared periodic cells."""

from itertools import product
import json
from pathlib import Path

import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle

CELLS = Path(__file__).with_name("structure_cells.json")


def read_cell(key):
    return json.loads(CELLS.read_text())["cells"][key]


def camera(azimuth=-32, elevation=22):
    phi, theta = np.deg2rad([azimuth, elevation])
    right = np.array([np.cos(phi), -np.sin(phi), 0])
    up = np.array([np.sin(phi)*np.sin(theta), np.cos(phi)*np.sin(theta), np.cos(theta)])
    return np.array([right, up, np.cross(right, up)])


def periodic_neighbours(cell):
    lattice, basis = np.array(cell["lattice"]), np.array(cell["fractional"])
    offsets = np.array(list(product(range(-2, 3), repeat=3)))
    image = (basis[None, :, :] + offsets[:, None, :]).reshape(-1, 3) @ lattice
    counts, distances = [], []
    for point in basis @ lattice:
        lengths = np.linalg.norm(image - point, axis=1)
        nearest = lengths[lengths > 1e-9].min()
        counts.append(int(np.sum(np.isclose(lengths, nearest, rtol=1e-7, atol=1e-9))))
        distances.append(float(nearest))
    return counts, distances


def geometry(key):
    cell = read_cell(key)
    lattice, basis = np.array(cell["lattice"]), np.array(cell["fractional"])
    if cell["prototype"] == "honeycomb":
        offsets = np.array([(i, j, 0) for i in range(4) for j in range(3)])
        points = (basis[None, :, :] + offsets[:, None, :]).reshape(-1, 3) @ lattice
        frame, edges = np.empty((0, 3)), []
    else:
        a = cell["conventional_a"]
        offsets = np.array(list(product(range(-2, 3), repeat=3)))
        points = (basis[None, :, :] + offsets[:, None, :]).reshape(-1, 3) @ lattice
        points = points[np.all((points >= -1e-9) & (points <= a+1e-9), axis=1)]
        points = np.unique(np.round(points, 12), axis=0)
        frame = np.array(list(product([0, a], repeat=3)))
        edges = [(i, j) for i in range(8) for j in range(i+1, 8)
                 if np.count_nonzero(frame[i] != frame[j]) == 1]
    counts, distances = periodic_neighbours(cell)
    expected = {"diamond": 4, "honeycomb": 3, "fcc": 12}[cell["prototype"]]
    assert counts == [expected] * len(basis), (key, counts)
    nearest = min(distances)
    bonds = [(i, j) for i in range(len(points)) for j in range(i+1, len(points))
             if np.isclose(np.linalg.norm(points[i]-points[j]), nearest, rtol=1e-7, atol=1e-9)]
    return cell, points, bonds, frame, edges, nearest


def sphere(ax, centre, radius, color, order):
    base = np.array(to_rgb(color))
    # Nested vector disks approximate directional surface lighting; no bitmap is embedded.
    for k, size in enumerate(np.linspace(1, .025, 36)):
        t = 1-size
        offset = radius * t * np.array([-.30, .38])
        shade = np.clip(base * (.46 + .52*t) + .46*t**2, 0, 1)
        ax.add_patch(Circle(centre[:2] + offset, radius*size, facecolor=shade,
                            edgecolor="none", zorder=order+k*.0001))


def draw_motif(fig, bounds_mm, key, color, azimuth=-32, elevation=22):
    cell, points, bonds, frame, edges, nearest = geometry(key)
    width, height = fig.get_size_inches()*25.4
    x, y, w, h = bounds_mm
    ax = fig.add_axes([x/width, y/height, w/width, h/height])
    ax.set_axis_off()
    view = camera(azimuth, elevation)
    centre = (points.min(axis=0)+points.max(axis=0))/2
    projected = (points-centre) @ view.T
    radius = nearest * {"fcc": .215, "honeycomb": .215, "diamond": .18}[cell["prototype"]]
    extent = np.max(np.abs(projected[:, :2]), axis=0)+1.1*radius
    scale = max(extent[0]/(w/h), extent[1])
    ax.set(xlim=(-scale*w/h, scale*w/h), ylim=(-scale, scale), aspect="equal")
    wire = (frame-centre) @ view.T if len(frame) else frame
    for i, j in edges:
        ax.plot(wire[[i, j], 0], wire[[i, j], 1], color="#8D969D", lw=.45,
                solid_capstyle="round", zorder=0)
    elements = [(float(point[2]), "atom", point) for point in projected]
    if cell["prototype"] != "fcc":
        for i, j in bonds:
            vector = projected[j]-projected[i]
            start = projected[i] + vector*radius/nearest
            end = projected[j] - vector*radius/nearest
            for t in range(8):
                segment = np.array([start+(end-start)*max(0, (t-.015)/8),
                                    start+(end-start)*min(1, (t+1.015)/8)])
                elements.append((float(segment[:, 2].mean()), "bond", segment))
    for order, (_, kind, value) in enumerate(sorted(elements, key=lambda item: item[0]), 1):
        if kind == "atom":
            sphere(ax, value, radius, color, order)
        else:
            for linewidth, shade, dz in ((1.45, "#6F7A82", 0), (.65, "#C5CBD0", .1)):
                ax.plot(value[:, 0], value[:, 1], lw=linewidth, color=shade,
                        solid_capstyle="butt", zorder=order+dz)
    return {"key": key, "element": cell["element"], "prototype": cell["prototype"],
            "display_atoms_including_boundary_images": len(points),
            "display_nearest_neighbour_pairs": len(bonds),
            "periodic_coordination": periodic_neighbours(cell)[0],
            "nearest_distance_bohr": nearest, "camera_degrees": [azimuth, elevation],
            "bounds_mm": list(bounds_mm), "vector_only": True,
            "size_interpretation": "periodic motif; not simulation-cell size"}
