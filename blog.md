---
layout: default
---
## Health economics with R

<ul style="margin:10px;padding:0">
  {% for post in site.posts %}
    <li class="blog-archive">
      <a href="{{ post.url }}">{{ post.title }}</a> <br>
       {{post.date | date: "%B %-d, %Y" }}
    </li>
  {% endfor %}
</ul>