[![Build Status](https://travis-ci.org/BodenmillerGroup/bbsnippets.svg?branch=master)](https://travis-ci.org/BodenmillerGroup/bbsnippets)

## Getting started

We are using Docsify [https://docsify.js.org] documentation site generator. We strongly recommend to put a `README.md`
file into the snippet's folder, as it'll be used during automatic documentation indexing. The first line in this file
is a title shown on a sidebar. Here is an example: [README.md](https://raw.githubusercontent.com/BodenmillerGroup/bbsnippets/master/snippets/imc/slide_overview/README.md)

All snippets are separated by domain-specific folders. Each such folder can contain multiple subfolders for different
tasks/topics. It is recommended to split examples in different languages (R, Python, etc.) into a separate directory,
so programming tools can manage these snippets accordingly.

## Inline source code

With docsify 4.6 it is now possible to embed any type of file. You can embed these files as video, audio, iframes,
or code blocks, and even Markdown files can even be embedded directly into the document.

For example, here embedded a Markdown file. You only need to do this:

```markdown
[filename](path/to/example.md ':include')
```

Normally, this will compiled into a link, but in docsify, if you add `:include` it will be embedded.

#### Embedded file type

Currently, file extension are automatically recognized and embedded in different ways.

This is a supported embedding type:

* **iframe** `.html`, `.htm`
* **markdown** `.markdown`, `.md`
* **audio** `.mp3`
* **video** `.mp4`, `.ogg`
* **code** other file extension

Of course, you can force the specified. For example, you want to Markdown file as code block embedded.

```markdown
[filename](path/to/example.md ':include :type=code')
```

##### Embedded code fragments
Sometimes you don't want to embed a whole file. Maybe because you need just a few lines but you want to compile and test
the file in CI.

```markdown
[filename](path/to/example.js ':include :type=code :fragment=demo')
```

In your code file you need to surround the fragment between `/// [demo]` lines (before and after the fragment).  
Alternatively you can use `### [demo]`.


##### Tag attribute

If you embed the file as `iframe`, `audio` and `video`, then you may need to set the attributes of these tags.
You can check [MDN](https://developer.mozilla.org/en-US/docs/Web/HTML/Element/iframe) for these attributes.

```markdown
[cinwell website](https://cinwell.com ':include :type=iframe width=100% height=400px')
```

##### The code block highlight

Embedding any type of source code file, you can specify the highlighted language or automatically identify.

```markdown
[](path/to/example.html ':include :type=code text')
```


## Code conventions

We suggest authors of submitted examples to follow commonly accepted source code styling practices. There are
.editorconfig and .prettierrc config files in the root folder that can be used to apply such rules. For code formatting
these tools provide a strongly opinionated styling rules with good editor integration:

- Python: Black (https://black.readthedocs.io)
- R: styler (https://styler.r-lib.org/)
- Other languages: Prettier (https://prettier.io/)
