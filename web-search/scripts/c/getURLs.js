#!/usr/bin/env node

/*
Extract all URLs from a web page.
Usage: ./getURLs.js <base_url>
*/

const readline = require('readline');
const {JSDOM} = require('jsdom');
const {URL} = require('url');

// 1. Read the base URL from the command-line argument using `process.argv`.
let baseURL = '';
// __start_solution__
baseURL = process.argv[2];
// __end_solution__

if (baseURL.endsWith('index.html')) {
  baseURL = baseURL.slice(0, baseURL.length - 'index.html'.length);
} else {
  baseURL += '/';
}

const rl = readline.createInterface({
  input: process.stdin,
});

// __start_solution__
let htmlInput = '';
// __end_solution__
rl.on('line', (line) => {
  // 2. Read HTML input from standard input (stdin) line by line using the `readline` module.
  // __start_solution__
  htmlInput += line + '\n';
  // __end_solution__
});

rl.on('close', () => {
  // 3. Parse HTML using jsdom
  // __start_solution__
  const dom = new JSDOM(htmlInput).window.document;
  // __end_solution__

  // 4. Find all URLs:
  //  - select all anchor (`<a>`) elements) with an `href` attribute using `querySelectorAll`.
  //  - extract the value of the `href` attribute for each anchor element.
  // __start_solution__
  const urls = Array.from(dom.querySelectorAll('a[href]')).
      map((element) => element.getAttribute('href'));
    // __end_solution__
    // 5. Print each absolute URL to the console, one per line.
    // __start_solution__
  urls.forEach((url) => {
    const absoluteURL = new URL(url, baseURL);
    console.log(absoluteURL.toString());
  });
  // __end_solution__
});


