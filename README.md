# WOLF PACK
The project aims to bundle our in-house development of Mathematica packages for structural analysis and optimization

** START DATE **: January-16-2018

** AUTHORS **: `Luis Bahamonde` (luisb@email.sc.edu), `Zafer Gurdal`

** STATUS **: `Development` phase for version [0.0.1](http://keepachangelog.com/en/1.0.0/#what)

Note:
- The email serves as a communication point with the maintainer of the project
- The authors are ordered as in research articles
- The version is linked to a changelog if you decide to have one

---
## Environment

If you strive for reproducible research, your computational environment must be reproducible.

This section briefly explains the environment to setup to get the _software project_ code running. The environmental components may come from a variety of sources. Use an `[environment.yml](./ENVIRONMENTS.yml)` file for third-party open source packages.

You may also be reusing `McNAIRPy` projects. On-premise projects used or serviced can be given as in the table below:

| from | for  | purpose |
| ---- |--------- | ----------- |
| [materials](./../materials/README.md) | |  uses materials database for anisotropic configuration definitions |
| | [cantilever](./../data/cantilever/README.md) | provides material distribution data file |

The from and for columns identify the upstream or downstream dependency. Note that you also list `McNAIRPy` projects that your _software project_ services, even if your _software project_ can run without them.

The intent of the table headers is to read each project as, for example, from materials the _software project_ uses some functionality and for cantilever the project provides this other functionality.

You can also, additionally, provide a stack figure of the components used.

Finally, you may have commercial off the shelf programs, like Abaqus. Briefly state them as
> **NOTE** Abaqus is needed for this project


If you want your development environment to be reproducible too, list the tool chain you used. The table below may be of help

| tool | description |
| ---- | ----------- |
| [Microsoft Teams](https://teams.microsoft.com/l/channel/19%3ab2b22aa0fe2e4ecf8a3f27a19f4ceb65%40thread.skype/utilities?groupId=a8e97c28-831c-4d26-b349-6a4c57d67c78&tenantId=4b2a4b19-d135-420e-8bb2-b1cd238998cc) | channel for pull requests, wiki |
| [McNAIRPy](///P:/README.md) | network drive that hosts our remote repositories for collaboration|
| [Git](https://git-scm.com/) | version control system |
| [changelog](http://keepachangelog.com/en/1.0.0/) | changelog for release version management |
| [Markdown Preview Plus](https://chrome.google.com/webstore/detail/markdown-preview-plus/febilkbfcbhebfnokafefeacimjdckgl?hl=en-US) | markdown renderer for changelog and readme files |
| [conda](https://conda.io/docs/user-guide/tasks/manage-pkgs.html) | package manager |
| [unittest](https://docs.python.org/2/library/unittest.html) | Python unit testing framework |
| [sphinx](http://www.sphinx-doc.org/en/stable/) | automatic documentation Python tool |

Note that Python tools used in your toolchain should be already in `[environment.yml](./ENVIRONMENTS.yml)`, but you can list them anyway

---
## Summary

This section provides a lower level (say 10000 feet) description of the architecture of your software program. Is it object-oriented? then describe the class diagram very briefly. Is it a procedure? what are the main steps of that procedure and what relevant data flow is present.

If you have a block diagram of your software system architecture, insert it as a figure:

```
![figure](./figs/mariaCompositesArchitecture.png)
```

If you want your software project to be imported as an alias show us
```
import my_software_project as bob
```

---
## References

Main references used to develop the project using AIAA format

1. [Software Documentation](./SPHINX.md)
