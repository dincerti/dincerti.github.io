---
layout: default
---

<ul style="margin:10px;padding:0">
  {% for post in site.posts %}
    <li class="blog-archive">
      <a href="{{ post.url }}">{{ post.title }}</a> <br>
       {{post.date | date: "%B %-d, %Y" }}
    </li>
  {% endfor %}
</ul>

<div class="footer">
   <p> Check out <a href="https://www.r-bloggers.com/" target="_blank"> R-bloggers </a> for other R content! </p>
   <p> <i class="fa fa-rss"></i>  <a href="feed.xml">RSS feed</a> </p>
</div>