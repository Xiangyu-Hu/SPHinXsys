"""SPHinXsys GDB pretty-printers (Eigen3-style).

This file is intentionally written in the same style as Eigen's official
printers: register a single lookup function via `obj.pretty_printers.append(...)`.

Target: SPH::ParticleVariables / BaseParticles::all_discrete_variables_
which is a std::tuple of std::vector<SPH::DiscreteVariable<...>*>.

To use it:
* Create a directory and put the file as well as an empty __init__.py in 
  that directory.
* Create a ~/.gdbinit file, that contains the following:
    python
    import sys
    sys.path.insert(0, '/path/to/sphinxsys/source/dir')
    from tools.gdb.printers import register_SPHinXsys_printers
    register_SPHinXsys_printers(None)
    end

With this installed you can do in GDB:
  (gdb) p all_discrete_variables_
in VSCode,
  (gdb) water_block.getBaseParticles().all_discrete_variables_
  (gdb) -exec p water_block.getBaseParticles().all_discrete_variables_
then you will get a grouped list of variable names.

Author: Hong Zhu, aided by GPT-5.2
"""

from __future__ import annotations

import re
import gdb


# When we are inside our own pretty-printer, calling gdb.default_visualizer(val)
# on the same value would pick our printer again and break child iteration.
# We use this flag to temporarily disable our lookup during internal traversal
# so libstdc++/Eigen printers can provide the tuple/vector children reliably.
_DISABLE_OUR_LOOKUP = False


_TYPE_LABELS = [
    "UnsignedInt",
    "int",
    "Real",
    "Vec2d (Vecd in 2D builds)",
    "Mat2d (Matd in 2D builds)",
    "Vec3d (Vecd in 3D builds)",
    "Mat3d (Matd in 3D builds)",
    "Vec6d",
    "Mat6d",
    "VecMatGrad2d",
    "VecMatGrad3d",
]


def _unquote(s: str) -> str:
    if len(s) >= 2 and ((s[0] == '"' and s[-1] == '"') or (s[0] == "'" and s[-1] == "'")):
        return s[1:-1]
    return s


def _try_default_children(val: gdb.Value):
    """Use already-loaded STL pretty-printers (libstdc++) when available."""
    global _DISABLE_OUR_LOOKUP
    try:
        prev = _DISABLE_OUR_LOOKUP
        _DISABLE_OUR_LOOKUP = True
        vis = gdb.default_visualizer(val)
    except Exception:
        vis = None
    finally:
        _DISABLE_OUR_LOOKUP = prev
    if vis is None:
        return None
    try:
        return list(vis.children())
    except Exception:
        return None


def _tuple_get_n_fallback(tuple_val: gdb.Value, n: int) -> gdb.Value:
    """Fallback tuple element access when no STL printer is available."""

    def descend(v: gdb.Value, idx: int) -> gdb.Value:
        if idx == 0:
            for f in v.type.fields():
                if f.name == "_M_head_impl":
                    return v[f]
            return v[v.type.fields()[0]]

        tail = None
        for f in v.type.fields():
            if f.name == "_M_tail":
                tail = v[f]
                break
        if tail is None:
            fields = v.type.fields()
            tail = v[fields[-1]]
        return descend(tail, idx - 1)

    return descend(tuple_val, n)


def _iter_tuple_elements(tuple_val: gdb.Value):
    children = _try_default_children(tuple_val)
    if children is not None:
        for _, elem in children:
            yield elem
        return

    i = 0
    while True:
        try:
            yield _tuple_get_n_fallback(tuple_val, i)
            i += 1
        except Exception:
            return


def _iter_vector_elements(vec_val: gdb.Value):
    children = _try_default_children(vec_val)
    if children is not None:
        for _, elem in children:
            yield elem
        return

    # Fallback for libstdc++ std::vector layout.
    try:
        start = vec_val["_M_impl"]["_M_start"]
        finish = vec_val["_M_impl"]["_M_finish"]
        it = start
        while it != finish:
            # element type is pointer, so dereference yields the pointer's target.
            yield it.dereference()
            it = it + 1
    except Exception:
        return


def _extract_entity_name(obj: gdb.Value):
    """Try to extract SPH::Entity::name_ from a (possibly derived) object."""
    # Prefer the public accessor; member name_ is protected and may not be
    # accessible depending on GDB settings / MI frontend.
    try:
        return _unquote(str(obj["Name"]()))
    except Exception:
        pass

    try:
        return _unquote(str(obj["name_"]))
    except Exception:
        pass

    # Search base classes (Entity may be a base).
    t = obj.type
    try:
        fields = t.fields()
    except Exception:
        fields = []

    for f in fields:
        if not getattr(f, "is_base_class", False):
            continue
        try:
            base_obj = obj.cast(f.type)
            return _unquote(str(base_obj["name_"]))
        except Exception:
            continue

    return None


class SPHParticleVariablesPrinter:
    """Pretty printer for SPH::ParticleVariables (all_discrete_variables_).

    This printer is intentionally *expandable* (tuple -> vectors -> variables)
    rather than returning a large single string.
    """

    def __init__(self, val: gdb.Value):
        self.val = val

    def to_string(self) -> str:
        return "std::tuple (SPH ParticleVariables / all_discrete_variables_)"

    def children(self):
        idx = 0
        for elem in _iter_tuple_elements(self.val):
            label = _TYPE_LABELS[idx] if idx < len(_TYPE_LABELS) else f"type_index_{idx}"
            yield (f"[{idx}] {label}", elem)
            idx += 1

    def display_hint(self):
        # VSCode renders this as an expandable object.
        return "tuple"


def _vector_size_capacity(vec_val: gdb.Value):
    """Return (size, capacity) for libstdc++ std::vector when possible."""
    try:
        start = vec_val["_M_impl"]["_M_start"]
        finish = vec_val["_M_impl"]["_M_finish"]
        end = vec_val["_M_impl"]["_M_end_of_storage"]
        size = int(finish - start)
        cap = int(end - start)
        return size, cap
    except Exception:
        # Best-effort fallback via STL visualizer length.
        try:
            children = _try_default_children(vec_val)
            if children is not None:
                return len(children), None
        except Exception:
            pass
    return None, None


def _extract_discrete_variable_data_size(obj: gdb.Value):
    try:
        return int(obj["data_size_"])
    except Exception:
        pass
    # Fallback to method call in case private member isn't accessible.
    try:
        return int(obj["getDataSize"]())
    except Exception:
        return None


class SPHDiscreteVariablePtrPrinter:
    """Pretty printer for SPH::DiscreteVariable<...>* pointers."""

    def __init__(self, val: gdb.Value):
        self.val = val

    def to_string(self) -> str:
        try:
            if int(self.val) == 0:
                return "SPH::DiscreteVariable* nullptr"
        except Exception:
            pass

        try:
            obj = self.val.dereference()
            name = _extract_entity_name(obj) or "<name_unavailable>"
            data_size = _extract_discrete_variable_data_size(obj)
            if data_size is None:
                return f"SPH::DiscreteVariable* {name}"
            return f"SPH::DiscreteVariable* {name} (data_size={data_size})"
        except Exception:
            return "SPH::DiscreteVariable* <unavailable>"

    def children(self):
        try:
            if int(self.val) == 0:
                return iter(())
        except Exception:
            return iter(())

        try:
            obj = self.val.dereference()
        except Exception:
            return iter(())

        # Expose name_ via Entity base when possible.
        try:
            yield ("name_", obj["name_"])
        except Exception:
            pass

        try:
            yield ("data_size_", obj["data_size_"])
        except Exception:
            pass

        # Keep raw pointer target visible as well.
        try:
            yield ("object", obj)
        except Exception:
            pass

    def display_hint(self):
        return "object"


class SPHDiscreteVariableVectorPrinter:
    """Pretty printer for std::vector<SPH::DiscreteVariable<...>*>."""

    def __init__(self, val: gdb.Value):
        self.val = val

    def to_string(self) -> str:
        size, cap = _vector_size_capacity(self.val)
        if size is None:
            return "std::vector<SPH::DiscreteVariable*>"
        if cap is None:
            return f"std::vector<SPH::DiscreteVariable*> size={size}"
        return f"std::vector<SPH::DiscreteVariable*> size={size} capacity={cap}"

    def children(self):
        # Render as an array so VSCode shows [0], [1], ...
        i = 0
        for p in _iter_vector_elements(self.val):
            # p is usually the pointer element value
            try:
                pv = p
                if pv.type.code != gdb.TYPE_CODE_PTR:
                    # In case the STL printer yielded a dereferenced value.
                    pv = pv.address
            except Exception:
                pv = p

            # Prefer showing variable name+data_size in the child key.
            key = str(i)
            try:
                if int(pv) != 0:
                    obj = pv.dereference()
                    name = _extract_entity_name(obj) or "<name_unavailable>"
                    data_size = _extract_discrete_variable_data_size(obj)
                    if data_size is None:
                        key = f"{i} ({name})"
                    else:
                        key = f"{i} ({name}, data_size={data_size})"
            except Exception:
                pass

            yield (key, pv)
            i += 1

    def display_hint(self):
        return "array"


# -------------------- Eigen-style dictionary lookup --------------------

pretty_printers_dict = {}


def build_sphinxsys_dictionary():
    # ParticleVariables is a tuple of vectors of DiscreteVariable pointers.
    pretty_printers_dict[
        re.compile(r"^std::tuple<.*SPH::DiscreteVariable.*>$")
    ] = lambda val: SPHParticleVariablesPrinter(val)


def _looks_like_particle_variables(type_str: str) -> bool:
    """Best-effort detection for BaseParticles::all_discrete_variables_.

    Different builds may show typedefs such as
    SPH::DataContainerAddressAssemble<SPH::DiscreteVariable,...> instead of
    the expanded std::tuple<...>.
    """
    s = type_str.replace(" ", "")
    if "DiscreteVariable" not in s:
        return False
    # The underlying structure is a tuple of vectors of DiscreteVariable pointers.
    if "std::tuple<" in s:
        return True
    if "DataContainerAddressAssemble<" in s:
        return True
    if "ParticleVariables" in s:
        return True
    return False


def lookup_function(val):
    """Look-up and return a pretty-printer that can print val."""
    if _DISABLE_OUR_LOOKUP:
        return None

    t = val.type
    if t.code == gdb.TYPE_CODE_REF:
        t = t.target()

    t = t.unqualified().strip_typedefs()

    # Prefer the full type string; tags are inconsistent across toolchains.
    typename = str(t)
    typename_nospace = typename.replace(" ", "")

    # Fast-path: ParticleVariables is a std::tuple of vectors.
    if _looks_like_particle_variables(typename):
        return SPHParticleVariablesPrinter(val)

    # Pointers to SPH::DiscreteVariable<...> only (avoid over-matching pointers
    # to vectors/tuples that merely contain DiscreteVariable in their template args).
    try:
        if t.code == gdb.TYPE_CODE_PTR:
            target_str = str(t.target().unqualified().strip_typedefs())
            target_nospace = target_str.replace(" ", "")
            if "SPH::DiscreteVariable<" in target_nospace and "std::vector<" not in target_nospace:
                return SPHDiscreteVariablePtrPrinter(val)
    except Exception:
        pass

    # std::vector<...DiscreteVariable<...>*...> (require the type itself is a vector)
    try:
        if typename_nospace.startswith("std::vector<") and "DiscreteVariable<" in typename_nospace:
            if "*>" in typename_nospace or "*const>" in typename_nospace or "*const,>" in typename_nospace:
                return SPHDiscreteVariableVectorPrinter(val)
            if "*" in typename_nospace:
                return SPHDiscreteVariableVectorPrinter(val)
    except Exception:
        pass

    for regex in pretty_printers_dict:
        if regex.search(typename):
            return pretty_printers_dict[regex](val)

    # Some libstdc++ builds include extra spaces; try a no-space match.
    typename_nospace = typename.replace(" ", "")
    for regex in pretty_printers_dict:
        if regex.search(typename_nospace):
            return pretty_printers_dict[regex](val)

    return None


def register_SPHinXsys_printers(obj):
    """Register SPHinXsys pretty-printers with objfile Obj (Eigen3 style)."""
    if obj is None:
        obj = gdb
    # Order matters: libstdc++ registers a generic std::tuple printer.
    # If we append, the tuple printer will match first and our printer will never run.
    # Insert at front to override tuple printing only for tuples that contain SPH::DiscreteVariable.
    pp = obj.pretty_printers
    # Remove existing occurrences (avoids duplicates, and ensures we can move to front).
    try:
        while True:
            pp.remove(lookup_function)
    except ValueError:
        pass
    except Exception:
        # If list ops behave oddly, fall back to best-effort insertion.
        pass

    try:
        pp.insert(0, lookup_function)
    except Exception:
        # As a last resort, append.
        pp.append(lookup_function)


build_sphinxsys_dictionary()