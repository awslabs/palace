// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// got this from the webpack docs
// https://webpack.js.org/configuration/output/
(function webpackUniversalModuleDefinition(root, factory) {
  // node
  if (typeof exports === "object" && typeof module === "object") {
    console.log("mod type 1");
    module.exports = factory(require("./glvis"));
  }
  // requirejs?
  else if (typeof define === "function" && define.amd) {
    console.log("mod type 2");
    define(["./glvis"], factory);
  }
  // unknown
  else if (typeof exports === "object") {
    console.warn("glvis: untested module load");
    exports["someLibName"] = factory(require("./glvis"));
  }
  // browser global
  else {
    root["glvis"] = factory(root["glvis"]);
  }
})(this, function (emglvis) {
  function rand_id() {
    const arr = new Uint8Array(10);
    window.crypto.getRandomValues(arr);
    return arr.reduce(
      (cur, next) => cur + next.toString(36).padStart(2, "0"),
      ""
    );
  }

  class State {
    constructor(div, width = 640, height = 480, canvas = undefined) {
      if (div === undefined) {
        throw "div cannot be undefined";
      }
      this.div_ = div;
      this.canvas_ = canvas;
      this.emglv_ = emglvis();
      this.emsetup_ = false;
      this.setupCanvas(width, height);
      this.new_stream_callbacks = [];
      // could also have an update_stream_callbacks
      // this.update_stream_callbacks = [];
    }

    setCanvasSize(width, height) {
      this.pixel_ratio_ = window.devicePixelRatio || 1;
      this.logical_width_ = width;
      this.logical_height_ = height;
      this.canvas_.style.width = `${width}px`;
      this.canvas_.style.height = `${height}px`;
      this.canvas_.width = Math.floor(width * this.pixel_ratio_);
      this.canvas_.height = Math.floor(height * this.pixel_ratio_);
      console.log(
        `dpr=${this.pixel_ratio_} new canvas sizes: ` +
          `style.width=${this.canvas_.style.width}, ` +
          `style.height=${this.canvas_.style.height}, ` +
          `width=${this.canvas_.width}, ` +
          `height=${this.canvas_.height}`
      );
    }

    async setSize(width, height) {
      this.setCanvasSize(width, height);
      var g = await this.emglv_;
      g.resizeWindow(width, height);
      g.sendExposeEvent();
    }

    setupCanvas(width, height) {
      if (this.canvas_ === undefined) {
        this.canvas_ = document.createElement("canvas");
        this.canvas_.id = `glvis-canvas-${rand_id()}`;
      }
      this.setCanvasSize(width, height);
      this.canvas_.innerHTML = "Your browser doesn't support HTML canvas";

      this.canvas_.oncontextmenu = function (e) {
        e.preventDefault();
        // prevent right-click from bubbling up to stuff like Jupyter Lab that
        // have a custom JS-based right-click menu we don't want opening when
        // someone is trying to zoom
        e.stopPropagation();
      };
      var that = this;
      this.canvas_.addEventListener("click", function () {
        that.div_.focus();
        return true;
      });

      this.div_.append(this.canvas_);
    }

    // only callable from resolved emglvis
    _setupEmGlvis(g) {
      if (this.emsetup_) {
        return;
      }
      g.setKeyboardListeningElementId(this.div_.id);
      g.setCanvasId(this.canvas_.id);
      g.canvas = this.canvas_;
    }

    // only callable from resolved emglvis
    _startVis(g) {
      if (this.emsetup_) {
        return;
      }
      this.emsetup_ = true;
      console.log("starting visualization loop");
      // needed again here... do we delete the SdlWindow somewhere in startVisualization?
      g.setCanvasId(this.canvas_.id);
      function iterVis(timestamp) {
        g.iterVisualization();
        window.requestAnimationFrame(iterVis);
      }
      window.requestAnimationFrame(iterVis);
    }

    async _display(data, is_update, ...args) {
      var g = await this.emglv_;
      if (Array.isArray(data)) {
        const f = is_update
          ? g.updateParallelStreams
          : g.displayParallelStreams;
        // TODO: this seems expensive, there must be a better way..
        var streams = new g.StringArray();
        data.forEach((s) => streams.push_back(s));
        const stat = f(streams, ...args);
        streams.delete();
        return stat;
      } else if (typeof data == "string") {
        const f = is_update ? g.updateStream : g.displayStream;
        return f(data, ...args);
      } else {
        throw `unsupported data type ${typeof data}`;
      }
    }

    async display(data) {
      var g = await this.emglv_;
      this._setupEmGlvis(g);

      await this._display(
        data,
        false,
        this.logical_width_,
        this.logical_height_
      );

      this.new_stream_callbacks.forEach((f) => f(this));
      this._startVis(g);
    }

    async update(data) {
      if (!this.emsetup_) {
        this.display(data);
        return;
      }
      if ((await this._display(data, true)) != 0) {
        console.log("unable to update stream, starting a new one");
        this.display(data);
      } else {
        console.log("updated stream");
      }
    }

    async sendKey(key, ctrl = false, shift = false, alt = false) {
      var key_code = undefined;
      if (typeof key === "string") {
        key_code = key.charCodeAt(0);
        console.log(`sending key '${key[0]}' (key_code=${key_code})`);
      } else if (typeof key === "number") {
        key_code = key;
        console.log(`sending key_code=${key_code}`);
      } else {
        throw "unsupported type";
      }
      var g = await this.emglv_;
      g.processKey(key_code, ctrl, shift, alt);
    }

    async sendKeyStr(keys) {
      var g = await this.emglv_;
      g.processKeys(keys);
    }

    async loadUrl(url) {
      var resp = await fetch(url);
      if (!resp.ok) {
        alert(`${url} doesn't exist`);
        return;
      }
      var text = await resp.text();
      if (text == "") {
        alert(`${url} has no content`);
        return;
      }
      console.log(`loading ${url}`);
      await this.display(text);
    }

    async loadStream(e) {
      const file = e.target.files[0];
      const data = await new Response(file).text();
      const extension = file.name.split('.').pop();
      if (["mesh", "vtk", "msh"].includes(extension)) {
        await this.display("mesh\n" + data);
      } else {
        await this.display(data);
      }
    }

    setTouchDevice(status) {
      // TODO TMS
      //this.emglv_.then((g) => { g.setTouchDevice(status); });
    }

    async getHelpString() {
      var g = await this.emglv_;
      return g.getHelpString();
    }

    // NOTE: the resulting data is just a view and the view will be invalid
    // after subsequent aclls to `getScreenBuffer` or `getPNGURL` since they
    // all use the same underlying buffer
    async getScreenBuffer(flip_y = false) {
      var g = await this.emglv_;
      return g.getScreenBuffer(flip_y);
    }

    async getScreenBufferAsCanvas() {
      let data = await this.getScreenBuffer(true);
      // idk why we need this but the `ImageData` ctor complains
      // when a Uint8Array is passed
      let clamped_data = new Uint8ClampedArray(data);
      const w = this.canvas_.width;
      const h = this.canvas_.height;
      let imdata = new ImageData(clamped_data, w, h);
      let can = document.createElement("canvas");
      can.width = w;
      can.height = h;
      let ctx = can.getContext("2d");
      ctx.putImageData(imdata, 0, 0);
      return can;
    }

    async getPNGAsB64() {
      let g = await this.emglv_;
      const data = g.getPNGByteArray();
      if (data === null) {
        throw "png data is null";
      }
      return btoa(
        data.reduce(function (data, byte) {
          return data + String.fromCharCode(byte);
        }, "")
      );
    }

    async getPNGURL() {
      return "data:image/png;base64," + (await this.getPNGAsB64());
    }

    async openScreenshotInTab() {
      const url = await this.getPNGURL();
      let tab = window.open("about:blank");
      tab.location.href = url;
    }

    async saveScreenshot(name = "glvis.png") {
      const url = await this.getPNGURL();
      let link = document.createElement("a");
      link.download = name;
      link.href = url;
      link.click();
    }

    // callbacks: f(State) -> void
    registerNewStreamCallback(f) {
      this.new_stream_callbacks.push(f);
    }
  }

  return {
    emglvis: emglvis,
    State: State,
    info: function () {
      console.log("hi from GLVis!");
    },
    rand_id: rand_id,
  };
});
