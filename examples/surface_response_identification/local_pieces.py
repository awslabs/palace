#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The LOCAL METAL PIECES of a translational cross-section (identification fix round DESIGN A.1 (c) / A.3 (b) / C.3;
decisions 1058 MAJOR-1, 1064 (4), 1096 (8)): ONE rule shared BY NAME with the C++ `LocalPiecesOfParallelEdges` /
`LocalPieceOfOffset` (palace/models/surfaceresponseoperator.cpp), used by the coupon generators (R0-b: one conductor
reference per local piece, Edges[].LocalPiece) and asserted equal to the C++ by U-PIECES (the C++ library reader
refuses a model whose Edges[].LocalPiece is not the piece its GapDirections bound).

Edges sorted by offset u; GapDirection +1 = the gap lies on the +u side of the edge (metal on its -u side), -1 = the
gap on the -u side (metal on +u). Edge k (gap -1) and edge k+1 (gap +1) bound ONE metal piece (the metal between
them); the metal on the -u side of the first edge (GapDirection +1) and on the +u side of the last edge (GapDirection
-1) are pieces reaching the box side. Piece ids 1..P in increasing u. A CPW cross-section (+1, -1, +1, -1) = 3 pieces:
ground-left / strip / ground-right; a gap (+1, -1) = 2; a strip (-1, +1) = 1."""


def local_pieces_of_parallel_edges(offsets, gap_directions):
    """-> [piece id per edge] (ints 1..P); ValueError (fail closed) on a non-cross-section, as the C++ MFEM_VERIFYs."""
    offsets = [float(u) for u in offsets]
    gap_directions = [int(g) for g in gap_directions]
    if not offsets or len(offsets) != len(gap_directions):
        raise ValueError("local_pieces_of_parallel_edges requires one GapDirection per edge offset")
    pieces = [0] * len(offsets)
    piece = 1 if gap_directions[0] == 1 else 0   # the piece left of the first edge exists iff its gap lies on its +u side
    for k in range(len(offsets)):
        if gap_directions[k] not in (1, -1):
            raise ValueError("local_pieces_of_parallel_edges: GapDirection must be -1 or 1")
        if k > 0 and not offsets[k] > offsets[k - 1]:
            raise ValueError("local_pieces_of_parallel_edges: edge offsets must be strictly increasing")
        if gap_directions[k] == 1:
            if k > 0 and gap_directions[k - 1] != -1:
                raise ValueError(f"local_pieces_of_parallel_edges: edges {k - 1} and {k} both have their gap on the +u side "
                                 "with no metal between them (not a cross-section)")
            pieces[k] = piece
        else:
            if k > 0 and gap_directions[k - 1] != 1:
                raise ValueError(f"local_pieces_of_parallel_edges: edges {k - 1} and {k} both have their gap on the -u side "
                                 "with no gap between them (not a cross-section)")
            piece += 1
            pieces[k] = piece
    return pieces


def local_piece_of_offset(offsets, gap_directions, pieces, u):
    """The piece whose metal interval holds u (0 when u lies in a gap): the C++ LocalPieceOfOffset."""
    offsets = [float(x) for x in offsets]
    if not (len(offsets) == len(gap_directions) == len(pieces)):
        raise ValueError("local_piece_of_offset requires the edges' offsets, gap directions and pieces")
    for k in range(len(offsets)):
        if int(gap_directions[k]) == 1:
            from_left = k == 0 or u >= offsets[k - 1]
            if from_left and u <= offsets[k]:
                return pieces[k]
        else:
            to_right = k + 1 == len(offsets) or u <= offsets[k + 1]
            if u >= offsets[k] and to_right:
                return pieces[k]
    return 0


def piece_count(pieces):
    return max(pieces) if pieces else 0


LOCAL_PIECES_MODEL_NAME_SUFFIX = "-pieces"


def local_pieces_model_name(base_name, pieces):
    """The piece-form model Name = <the one-reference form's Name>-pieces<P> (decision 1101; review d1099 MINOR-3):
    a library holding both forms of one key (DESIGN A.3 (e): the old entries kept for the ES record form) carries two
    distinct Names. The generator applies it to its default Name; the planner passes the suffixed --model-name."""
    return f"{base_name}{LOCAL_PIECES_MODEL_NAME_SUFFIX}{piece_count(pieces)}"


def metal_intervals_by_piece(offsets, gap_directions, pieces, radius):
    """{piece: (left, right)} of every piece's metal interval, the box-side pieces reaching offsets[0] - radius /
    offsets[-1] + radius (the coupon's half-width beyond the outer edges)."""
    offsets = [float(u) for u in offsets]
    bounds = [offsets[0] - float(radius)] + offsets + [offsets[-1] + float(radius)]
    out = {}
    for k in range(len(offsets)):
        if int(gap_directions[k]) == 1:
            left = bounds[k]          # the previous edge (or the box side)
            right = offsets[k]
        else:
            left = offsets[k]
            right = bounds[k + 2]     # the next edge (or the box side)
        p = pieces[k]
        lo, hi = out.get(p, (left, right))
        out[p] = (min(lo, left), max(hi, right))
    return out


def mesh_label(edge):
    """The label the coupon MESH and configs carry for a library entry's edge (Edges[] / ContextEdges[]): its LocalPiece
    in the piece form (R0-b: the mesher tags the metal per piece), its Conductor (the key label) otherwise. The consumers
    that compare the entry with the mesher's CSVs / the coupon configs read this, never Conductor alone."""
    return int(edge["LocalPiece"]) if edge.get("LocalPiece") is not None else int(edge["Conductor"])
