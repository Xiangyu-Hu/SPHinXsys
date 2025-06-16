---
layout: post
title:  "Debug SPHinXsys in VS Code with vcpkg installed library"
date:   2025-04-23
categories: sphinxsys update
---

You can do line-by-line debug of SPHinXsys in VS Code even into the vcpkg installed library.

By default, the vcpkg installed debug library is in the directory of `~\vcpkg\installed\x64-windows\debug\lib` 
or `~\vcpkg\installed\x64-linux\debug\lib`.
The source code is not included in the library,
but you can find it in the directory of `~\vcpkg\buildtrees\<library_name>\src\<version>`.
However, vcpkg usually does not keep the source code after the library is built.
To do line-by-line debug, you need to

1. Remove the library fully by `./vcpkg remove <library_name> --recurse`.

2. Reinstall the library by `./vcpkg install <library_name> --editable`.

After that, GDB will find the source files and the debug can be carried on.

<script src="https://giscus.app/client.js"
        data-repo="Xiangyu-Hu/SPHinXsys"
        data-repo-id="MDEwOlJlcG9zaXRvcnkxODkwNzAxNDA="
        data-category="Announcements"
        data-category-id="DIC_kwDOC0T7PM4CPNAR"
        data-mapping="pathname"
        data-strict="0"
        data-reactions-enabled="1"
        data-emit-metadata="0"
        data-input-position="bottom"
        data-theme="light"
        data-lang="en"
        crossorigin="anonymous"
        async>
</script>
